# 加载包
library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(stringr)
library(DoubletFinder)

# =====================【只在这里修改路径！】=====================
path_10x_data  <- "./02seurat/01data"
path_save_rdata<- ".data/data"
path_plot_out  <- "./01pre_deal_data"
path_allcell   <- "./00_all"
# ==============================================================

# Filter Sample Folders
folders <- list.files(path = path_10x_data, pattern = '[0124]$')
folders
if(length(folders) ==0) stop("ERROR: No folder matching [0124]$ was found under path_10x_data; please check the path!")

# Batch create Seurat objects, min.cells=3, min.features=200
sceList <- lapply(folders, function(folder_name){
  CreateSeuratObject(
    counts = Read10X(data.dir = file.path(path_10x_data, folder_name)),
    project = folder_name,
    min.cells = 3,
    min.features = 200
  )
})

seurat_filterDouble <- list()

# DoubletFinder cycle
for(n in seq_along(sceList)){
  obj <- sceList[[n]]
  # Mitochondrial proportion
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  # Basic quality control filtering
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 30)
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method='vst')
  obj <- obj %>% ScaleData(verbose = FALSE) %>% RunPCA(npcs = 40, verbose = FALSE)
  
  # Doublet parameter search
  sweep.res.list_kidney <- paramSweep_v3(obj, PCs = 1:40)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  mpK <- as.numeric(as.vector(bcmvn_kidney$pK[which.max(bcmvn_kidney$BCmetric)]))
  
  annotations <- obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(0.075 * ncol(obj))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  obj_doub <- doubletFinder_v3(obj, PCs = 1:40, pN = 0.25, pK = mpK, 
                               nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  seurat_filterDouble[[n]] <- obj_doub
}

# Save the list of doublet results.
save(seurat_filterDouble, file = file.path(path_save_rdata, "01seurat_filterDouble.Rdata"))

# =========【Key Optimization: Automatically extract DataFrame columns—no more manually writing out 33 lines!】=========
# Iteratively remove doublets, retaining only singlets.
seurat_singlet <- lapply(seurat_filterDouble, function(obj){
  # Automatically identify meta columns starting with 'DF.classifications'
  df_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)
  if(length(df_col) !=1) stop("找不到DF列或者多个DF列！")
  obj_sub <- subset(obj, subset = .data[[df_col]] == "Singlet")
  return(obj_sub)
})

# Merge all samples
sce.big <- merge(
  x = seurat_singlet[[1]],
  y = seurat_singlet[-1],
  add.cell.ids = folders,
  project = "CRC"
)

# ================= Grouping mapping table; replaces excessively long if-else chains; modify groupings only here ================
# 1. TWO_Curativeeffect: Drug efficacy grouping
group_effect <- c(
  C100 = "Resistant", C102 = "Resistant", C104 = "Resistant",
  C110 = "Sensitive", C12  = "Resistant", C120 = "Unknown",
  C122 = "Sensitive", C130 = "Sensitive", C132 = "Sensitive",
  C134 = "Sensitive", C20  = "Sensitive", C30  = "Resistant",
  C32  = "Resistant", C40  = "Resistant", C42  = "Resistant",
  C44  = "Resistant", C50  = "Resistant", C52  = "Resistant",
  C54  = "Resistant", C60  = "Resistant", C62  = "Resistant",
  C64  = "Resistant", C70  = "Sensitive", C80  = "Resistant",
  C90  = "Sensitive", C92  = "Sensitive", C94  = "Sensitive",
  N114 = "Sensitive", N134 = "Sensitive", N44  = "Resistant",
  N64  = "Resistant"
)

# 2. TimePoints (Time/Sample Type Grouping)
group_time <- c(
  C20 = "CRC_0", C30 = "CRC_0", C40 = "CRC_0", C50 = "CRC_0",
  C60 = "CRC_0", C70 = "CRC_0", C80 = "CRC_0", C90 = "CRC_0",
  C100= "CRC_0", C110= "CRC_0", C120= "CRC_0", C130= "CRC_0",
  C32 = "CRC_2", C12 = "CRC_2", C42 = "CRC_2", C52 = "CRC_2",
  C62 = "CRC_2", C92 = "CRC_2", C102= "CRC_2", C122= "CRC_2",
  C132= "CRC_2", C44 = "CRC_4", C54 = "CRC_4", C64 = "CRC_4",
  C94 = "CRC_4", C114= "CRC_4", C104= "CRC_4", C134= "CRC_4",
  C34 = "CRC_4", N44 = "CRC_Normal", N114="CRC_Normal",
  N134="CRC_Normal", N64 = "CRC_Normal"
)
# ==========================================================================

# Extract sample ID (the first part of the cell name before the underscore)
sample_id <- str_split(colnames(sce.big), pattern = "\\_", simplify = TRUE)[,1]

# Matching group; defaults to Sensitive / CRC_Normal if not in lookup.
TWO_Curativeeffect <- ifelse(is.na(group_effect[sample_id]), "Sensitive", group_effect[sample_id])
TimePoints         <- ifelse(is.na(group_time[sample_id]), "CRC_Normal", group_time[sample_id])

sce.big <- AddMetaData(sce.big, metadata = TWO_Curativeeffect, col.name = 'TWO_Curativeeffect')
sce.big <- AddMetaData(sce.big, metadata = TimePoints, col.name = 'TimePoints')

# View grouped statistics
table(sce.big@meta.data$TimePoints, sce.big@meta.data$TWO_Curativeeffect)

# High-Quality Violin Plot
pdf(file.path(path_plot_out, "01sub_single_information.pdf"))
print(VlnPlot(sce.big, features =c("nFeature_RNA","nCount_RNA","percent.mt"), ncol =3))
dev.off()

# Integrated Dimensionality Reduction and Clustering
sce.big1 <- FindVariableFeatures(sce.big, selection.method='vst', nfeatures = 5000)
sce.big1 <- sce.big1 %>% ScaleData(verbose = FALSE) %>% RunPCA(npcs = 40, verbose = FALSE)
sce.big1 <- RunUMAP(sce.big1, reduction = "pca", dims = 1:35)
sce.big1 <- FindNeighbors(sce.big1, reduction = "pca", dims = 1:35)
sce.big1 <- FindClusters(sce.big1, resolution = 0.8)
sce.big1

# Save the merged RDS file
saveRDS(sce.big1, file = file.path(path_save_rdata, "02CRC_ALL_single.rds"))

# Load the Harmony object and extract the intersecting cells.
allcell <- readRDS(file.path(path_allcell, "all_final_CRC.rds"))
common_cells <- intersect(colnames(allcell), rownames(sce.big1@meta.data))
allcell_sub <- subset(allcell, cells = common_cells)
saveRDS(allcell_sub, file = file.path(path_allcell, "03CRC_ALL_single_harmony.rds"))