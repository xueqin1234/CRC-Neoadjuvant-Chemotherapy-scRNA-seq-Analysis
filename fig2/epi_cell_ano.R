#=========================================================================
# Script: Epithelial cell subclustering and cell proportion enrichment analysis
# Description: Load CRC seurat object, annotate epithelial subtypes,
#              perform UMAP clustering, calculate cell proportion and CMH test
# Working directory: /hanlab/xueqinliu/test
#=========================================================================
library(Seurat)
library(harmony)
library(cowplot)
library(ggplot2)
library(sctransform)
library(future)
library(glmGamPoi)
library(stringr)
library(dplyr)
library(DoubletFinder)
library(vcdExtra)
library(vcd)
library(pheatmap)

#==================== Global setting & Input path ====================
# ========= 全部工作目录切换为 /hanlab/xueqinliu/test =========
root_dir <- "./test"
raw_data_dir <- "./test/data"
external_rds_dir <- "./test/00data"
# Load main CRC Seurat object
allcell <- readRDS(file.path(raw_data_dir, "03CRC_ALL_single_harmony.rds"))
allcell@meta.data$cell_name <- rownames(allcell@meta.data)

# Extract sample ID from cell barcode
sample_id <- sapply(colnames(allcell), function(x) unlist(strsplit(x, "\\_"))[1])

#==================== Add metadata: Group_Curativeeffect (3 response groups) ====================
# Three-category treatment response: Resistant / Sensitive / Intermediate / Unknown
group_3 <- case_when(
  sample_id %in% c('C100','C102','C104','C12','C50','C52','C54','C60','C62','C64','N64') ~ "Resistant",
  sample_id %in% c('C110','C122','C130','C132','C134','C20','C70','C90','C92','C94','N114','N134') ~ "Sensitive",
  sample_id %in% c('C30','C32','C40','C42','C44','C80','N44') ~ "Intermediate",
  sample_id == 'C120' ~ "Unknown",
  TRUE ~ "Sensitive"
)
allcell <- AddMetaData(allcell, metadata = group_3, col.name = 'Group_Curativeeffect')

#==================== Add metadata: Curativeeffect (binary response) ====================
# Binary treatment response grouping
group_2 <- case_when(
  sample_id %in% c('C50','C60','C100','C62','C102','C54','C64','C104','C52','C80','N64') ~ "Resistant",
  sample_id == 'C120' ~ "Unknown",
  TRUE ~ "Sensitive"
)
allcell <- AddMetaData(allcell, metadata = group_2, col.name = 'Curativeeffect')

#==================== Rename cell type identities & filter cells ====================
allcell1 <- RenameIdents(
  allcell,
  `0` = "T cells",
  `1` = "Epithelial cells",
  `2` = "T cells",
  `3` = "B cells",
  `4` = "Epithelial cells",
  `5` = "Myeloid",
  `6` = "Stromal cells",
  `7` = "T cells",
  `8` = "B cells",
  `9` = "B cells",
  `10` = "Epithelial cells",
  `11` = "Myeloid",
  `12` = "Stromal cells",
  `13` = "Epithelial cells",
  `14` = "Basal",
  `15` = "Stromal cells",
  `16` = "Mast",
  `17` = "Epithelial cells",
  `18` = "Myeloid",
  `19` = "Basal",
  `20` = "Epithelial cells",
  `21` = "B cells",
  `22` = "Basal",
  `23` = "Myeloid",
  `24` = "Epithelial cells",
  `25` = "Epithelial cells",
  `26` = "Epithelial cells",
  `27` = "T cells",
  `28` = "Stromal cells",
  `29` = "B cells",
  `30` = "Acinar cells",
  `31` = "Stromal cells",
  `32` = "Myeloid",
  `33` = "Epithelial cells",
  `34` = "T cells",
  `35` = "Epithelial cells",
  `36` = "T cells",
  `37` = "Stromal cells",
  `38` = "B cells",
  `39` = "B cells",
  `40` = "Myeloid"
)
table(Idents(allcell1))

# Remove Basal, Acinar cells and sample C114
allcell2 <- subset(allcell1, idents = c("Basal", "Acinar cells"), invert = TRUE)
allcell2 <- subset(allcell2, orig.ident != "C114")

# Set factor order for major cell types
Idents(allcell2) <- factor(
  Idents(allcell2),
  levels = c("Epithelial cells","Stromal cells","T cells","B cells","Myeloid","Mast")
)
table(Idents(allcell2))

# Subset only epithelial cells for downstream analysis
Epi <- subset(allcell2, idents = "Epithelial cells")

#==================== Add CMS subtype annotation ====================
# 依赖文件 cell_CMS.rds 放到 /hanlab/xueqinliu/test
cell_CMS <- readRDS("cell_CMS.rds")
Epi <- subset(Epi, cells = cell_CMS$cell_ID)

cell_CMS$cell_name <- cell_CMS$cell_ID
input_meta <- left_join(Epi@meta.data, cell_CMS, by = "cell_name")
Epi@meta.data$CMS <- input_meta$CMS

#==================== Add malignant / normal / non-malignant epithelial label (tomur) ====================
# 依赖文件 all_malignant_normal_epi.rds 放到 /hanlab/xueqinliu/test
malignant_ref <- readRDS("all_malignant_normal_epi.rds")
malignant_ref <- subset(malignant_ref, tomur == "malignant")
malignant_cell <- rownames(malignant_ref@meta.data)

# Subset malignant epithelial cells from Epi object
malignant_epi <- subset(Epi, cells = malignant_cell)
malignant_meta <- malignant_epi@meta.data
malignant_meta$tomur <- "malignant"

# Normal epithelium (CRC_Normal timepoint)
normal_epi <- subset(Epi, TimePoints == "CRC_Normal")
normal_epi_meta <- normal_epi@meta.data
normal_epi_meta$tomur <- "normal_epi"

# Non-malignant epithelium (exclude malignant and CRC_Normal)
nomalignant <- subset(Epi, cells = malignant_cell, invert = TRUE)
nomalignant <- subset(nomalignant, TimePoints != "CRC_Normal")
nomalignant_meta <- nomalignant@meta.data
nomalignant_meta$tomur <- "nomalignant"

# Merge tomur metadata back into Epi
meta_tomur <- bind_rows(malignant_meta, normal_epi_meta, nomalignant_meta)
input_meta <- left_join(Epi@meta.data, meta_tomur, by = "cell_name")
Epi@meta.data$tomur <- input_meta$tomur

#==================== Recreate Seurat object and reprocess epithelium ====================
Epi <- CreateSeuratObject(counts = Epi@assays$RNA@counts, meta.data = Epi@meta.data)

# Remove mitochondrial genes
Epi[["percent.mt"]] <- PercentageFeatureSet(Epi, pattern = "^MT-")
mt.index <- grep("^MT-", rownames(Epi), value = FALSE)
Epi <- Epi[-mt.index, ]

# Standard Seurat workflow: normalization -> variable genes -> PCA -> UMAP -> clustering
Epi <- NormalizeData(Epi, verbose = FALSE)
Epi <- FindVariableFeatures(Epi, selection.method = 'vst', nfeatures = 2000)
var_gene <- VariableFeatures(Epi)
Epi <- ScaleData(Epi, features = var_gene)
Epi <- RunPCA(Epi, features = var_gene, npcs = 80, verbose = FALSE)

Epi <- RunUMAP(Epi, reduction = "pca", dims = 1:40, verbose = FALSE)
Epi <- FindNeighbors(Epi, reduction = "pca", dims = 1:40, verbose = FALSE)
Epi <- FindClusters(Epi, resolution = 0.6, verbose = FALSE)

# Plot UMAP figures，输出到test目录
pdf("01_2000_06_all_epi.pdf")
DimPlot(Epi, group.by = "TimePoints")
DimPlot(Epi, group.by = "Group_Curativeeffect")
DimPlot(Epi, group.by = "tomur")
DimPlot(Epi, reduction = "umap", group.by = "orig.ident", pt.size = .1)
DimPlot(Epi, reduction = "umap", label = TRUE, pt.size = .1)
dev.off()

# Find all marker genes for epithelial clusters，csv输出到test
df_marker <- FindAllMarkers(Epi, verbose = FALSE)
write.csv(df_marker, file = "01_2000_06_allepi_maker.csv", row.names = FALSE)

#==================== Rename epithelial subclusters & filter non-epithelial contaminations ====================
Epit <- RenameIdents(
  Epi,
  `0` = "Mature enterocyte type II",
  `1` = "Stem/TA-like",
  `2` = "Goblet",
  `3` = "Unknown",
  `4` = "Stem/Progenitor cells",
  `5` = "Immature Goblet",
  `6` = "Immature Goblet",
  `7` = "Mature enterocyte type II",
  `8` = "Stem/TA-like",
  `9` = "Stem/Progenitor cells",
  `10` = "Paneth cells",
  `11` = "Stem/Progenitor cells",
  `12` = "Stem/Progenitor cells",
  `13` = "Mature enterocyte type II",
  `14` = "Mature enterocyte type II",
  `15` = "Stem/Progenitor cells",
  `16` = "Mature enterocyte type I",
  `17` = "Stem/TA-like",
  `18` = "Mature enterocyte type II",
  `19` = "Paneth cells",
  `20` = "Stem/Progenitor cells",
  `21` = "Stem/Progenitor cells",
  `22` = "Stem/Progenitor cells",
  `23` = "Enteroendocrine",
  `24` = "Stem/TA-like",
  `25` = "Mature enterocyte type II",
  `26` = "Tuft",
  `27` = "Stem/TA-like",
  `28` = "B",
  `29` = "Immature Goblet",
  `30` = "Unknown",
  `31` = "Stem/Progenitor cells",
  `32` = "Goblet",
  `33` = "Stromal cells"
)

# Remove contaminating stromal and B cells
Epit <- subset(Epit, idents = c("Stromal cells","B"), invert = TRUE)

# Set custom factor order for epithelial subtypes
Idents(Epit) <- factor(
  Idents(Epit),
  levels = c("Mature enterocyte type I","Mature enterocyte type II","Immature Goblet","Goblet",
             "Stem/TA-like","Stem/Progenitor cells","Paneth cells","Enteroendocrine","Tuft","Unknown")
)
table(Idents(Epit))

# UMAP plot with custom color，pdf输出test
pdf("01-1all_epi_cluster.pdf")
DimPlot(Epit, reduction = "umap", pt.size = 0.5,
        cols = c("#74b4f3","#77d1fc","#f4cee1","#f6aba6","#9fd1a8","#f9deb5","#a3d1e1","#7cb87a","#fad7ce","#d592be"))
dev.off()

Epit@meta.data$subcluster <- Idents(Epit)

#==================== Re-attach treatment response metadata to Epit object ====================
sample_id_epit <- sapply(colnames(Epit), function(x) unlist(strsplit(x, "\\_"))[1])

# 3-category response
group3_epit <- case_when(
  sample_id_epit %in% c('C100','C102','C104','C12','C50','C52','C54','C60','C62','C64','N64') ~ "Resistant",
  sample_id_epit %in% c('C110','C122','C130','C132','C134','C20','C70','C90','C92','C94','N114','N134') ~ "Sensitive",
  sample_id_epit %in% c('C30','C32','C40','C42','C44','C80','N44') ~ "Intermediate",
  sample_id_epit == 'C120' ~ "Unknown",
  TRUE ~ "Sensitive"
)
Epit <- AddMetaData(Epit, metadata = group3_epit, col.name = 'Group_Curativeeffect')

# Binary response
group2_epit <- case_when(
  sample_id_epit %in% c('C50','C60','C100','C62','C102','C54','C64','C104','C52','C80','N64') ~ "Resistant",
  sample_id_epit == 'C120' ~ "Unknown",
  TRUE ~ "Sensitive"
)
Epit <- AddMetaData(Epit, metadata = group2_epit, col.name = 'Curativeeffect')

#==================== Cell proportion statistics & CMH test ====================
# 不再单独setwd，直接在test目录输出
# Subset cells by timepoint, remove Unknown group
Epit_CRC_0 <- subset(Epit, TimePoints=="CRC_0" & Group_Curativeeffect!="Unknown")
Epit_CRC_2 <- subset(Epit, TimePoints=="CRC_2" & Group_Curativeeffect!="Unknown")
Epit_CRC_4 <- subset(Epit, TimePoints=="CRC_4" & Group_Curativeeffect!="Unknown")
Epit_CRC_Normal <- subset(Epit, TimePoints=="CRC_Normal" & Group_Curativeeffect!="Unknown")

# Count cell number matrix
specificity_CRC_0 <- table(Epit_CRC_0@meta.data$subcluster, Epit_CRC_0@meta.data$Group_Curativeeffect)
colnames(specificity_CRC_0) <- c("Intermediate_0","Resistant_0","Sensitive_0")

specificity_CRC_2 <- table(Epit_CRC_2@meta.data$subcluster, Epit_CRC_2@meta.data$Group_Curativeeffect)
colnames(specificity_CRC_2) <- c("Intermediate_2","Resistant_2","Sensitive_2")

specificity_CRC_4 <- table(Epit_CRC_4@meta.data$subcluster, Epit_CRC_4@meta.data$Group_Curativeeffect)
colnames(specificity_CRC_4) <- c("Intermediate_4","Resistant_4","Sensitive_4")

specificity_CRC_Normal <- table(Epit_CRC_Normal@meta.data$subcluster, Epit_CRC_Normal@meta.data$TimePoints)

# Combine count matrix
specificity <- cbind(specificity_CRC_0, specificity_CRC_2, specificity_CRC_4, specificity_CRC_Normal)

# CMH test and calculate observed/expected ratio
CMH_result <- CMHtest(specificity)
expected_mat <- chisq.test(specificity)$expected
ratio_obs_exp <- specificity / expected_mat

# Z-score normalization
df_t <- t(ratio_obs_exp)
z_mat <- scale(df_t, center = TRUE, scale = TRUE)
z_mat <- t(z_mat)

write.csv(z_mat, file = "Epit_all_R_obs_exp.csv")

# Plot pheatmap of z-score
annotation_col = data.frame(
  TimePoints = factor(rep(c("CRC_0","CRC_2","CRC_4","Normal"), c(3, 3, 3,1))),
  Group = factor(c(rep(c("Intermediate","Resistant","Sensitive"),3), rep("Normal",1)))
)
rownames(annotation_col) = colnames(z_mat)

ann_colors = list(
  TimePoints = c(CRC_0 = "#008000", CRC_2 = "#0000FF", CRC_4 = "#DC143C","Normal"="#BC8F8F"),
  Group = c(Intermediate="#2F4F4F",Resistant="#B8860B",Sensitive="#008B8B",Normal="#BC8F8F")
)

pdf("02Epit_Group_Curativeeffect_R_obs_exp.pdf", width=4, height=4)
pheatmap(z_mat, color = colorRampPalette(c("#0000cc","#f5f5f5","#cc6600","#cc0000"))(50),
         cluster_row = FALSE, cluster_col = FALSE,
         annotation_col = annotation_col, annotation_colors = ann_colors)
pheatmap(z_mat, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_row = FALSE, cluster_col = FALSE,
         annotation_col = annotation_col, annotation_colors = ann_colors)
dev.off()

# Save final epithelial object，rds保存到test目录
saveRDS(Epit, file.path(external_rds_dir, "all_epi_umap.rds"))
