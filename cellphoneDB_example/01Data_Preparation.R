#==================== Load packages ====================
library(Seurat)
library(dplyr)

#==================== Path Configuration (only modify root_dir here) ====================
root_dir <- "/hanlab/xueqinliu/test"
# Raw data path (original rds not modified)
raw_data_dir <- "/hanlab/xueqinliu/test/data"
external_rds_dir <- "/hanlab/xueqinliu/test/00data"
cpdb_root <- file.path(root_dir, "cellphoneDB_to_CRC/cellphoneDB_laoban/Intermediate")

# Auto create cellphoneDB subfolders
cpdb_crc0 <- file.path(cpdb_root, "CRC0")
cpdb_crc2 <- file.path(cpdb_root, "CRC2")
cpdb_crc4 <- file.path(cpdb_root, "CRC4")
dir.create(cpdb_crc0, recursive = TRUE, showWarnings = FALSE)
dir.create(cpdb_crc2, recursive = TRUE, showWarnings = FALSE)
dir.create(cpdb_crc4, recursive = TRUE, showWarnings = FALSE)

#==================== Load main CRC seurat object ====================
allcell <- readRDS(file.path(raw_data_dir, "03CRC_ALL_single_harmony.rds"))
allcell@meta.data$cell_name <- rownames(allcell@meta.data)

# Extract sample ID from cell barcode
sample_id <- sapply(colnames(allcell), function(x) unlist(strsplit(x, "\\_"))[1])

#==================== Add 3-category treatment response metadata: Group_Curativeeffect ====================
response_3 <- case_when(
  sample_id %in% c('C100','C102','C104','C12','C50','C52','C54','C60','C62','C64','N64') ~ "Resistant",
  sample_id %in% c('C110','C122','C130','C132','C134','C20','C70','C90','C92','C94','N114','N134') ~ "Sensitive",
  sample_id %in% c('C30','C32','C40','C42','C44','C80','N44') ~ "Intermediate",
  sample_id == 'C120' ~ "Unknown",
  TRUE ~ "Sensitive"
)
allcell <- AddMetaData(allcell, metadata = response_3, col.name = 'Group_Curativeeffect')

#==================== Add binary response metadata: Curativeeffect ====================
# NOTE: Original code sets almost all samples to Sensitive, keep original rule unchanged
response_2 <- case_when(
  sample_id %in% c() ~ "Resistant",
  sample_id == 'C120' ~ "Unknown",
  TRUE ~ "Sensitive"
)
allcell <- AddMetaData(allcell, metadata = response_2, col.name = 'Curativeeffect')

#==================== Rename cell identities & filter cells ====================
allcell1 <- RenameIdents(allcell, 
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
                         `40` = "Myeloid")
table(Idents(allcell1))

# Remove Basal, Acinar cells and sample C114
allcell2 <- subset(allcell1, idents = c("Basal", "Acinar cells"), invert = TRUE)
allcell2 <- subset(allcell2, orig.ident == "C114", invert = TRUE)

# Set factor order for cell types
Idents(allcell2) <- factor(Idents(allcell2), 
                           levels = c("Epithelial cells","Stromal cells","T cells","B cells","Myeloid","Mast"))
table(Idents(allcell2))

#==================== Load pre-annotated subcluster RDS files ====================
stromal <- readRDS(file.path(external_rds_dir, "stromal_umap_new.rds"))
tcell <- readRDS(file.path(external_rds_dir, "tcell_umap.rds"))
Bcell <- readRDS(file.path(external_rds_dir, "B_umap.rds"))
epi <- readRDS(file.path(external_rds_dir, "all_epi_nomalignant_malignant.rds"))
Myeloi <- readRDS(file.path(external_rds_dir, "01Myeloi_umap.rds"))
epi1 <- readRDS(file.path(external_rds_dir, "01_rank8_malignant.nmf.rds"))

epi1@meta.data$NMF_cluster <- paste("NMF", epi1$cluster, sep="_")

# Split epithelial cells into res / other based on NMF cluster
meta1 <- epi1@meta.data %>% 
  filter(NMF_cluster %in% c("NMF_1","NMF_4","NMF_5")) %>%
  select(cell_name, NMF_cluster, tomur)
meta1$subcluster <- "res"
meta1 <- meta1 %>% select(cell_name, subcluster, tomur)

meta2 <- epi1@meta.data %>% 
  filter(NMF_cluster %in% c("NMF_2","NMF_3","NMF_6","NMF_7","NMF_8")) %>%
  select(cell_name, NMF_cluster, tomur)
meta2$subcluster <- "other"
meta2 <- meta2 %>% select(cell_name, subcluster, tomur)

nmf_meta <- rbind(meta1, meta2)
colnames(nmf_meta) <- c("cell_name","subcluster","tumor")

# Non-malignant epithelial metadata
nomalignant <- subset(epi, tomur=="malignant", invert=TRUE)
nomalignant_meta <- nomalignant@meta.data %>% select(cell_name, subcluster, tomur)
colnames(nomalignant_meta) <- c("cell_name","subcluster","tumor")

# Stromal meta
stromal@meta.data$cell_name <- rownames(stromal@meta.data)
stromal_meta <- stromal@meta.data %>% select(cell_name, subcluster)
stromal_meta$tumor <- stromal_meta$subcluster

# T cell meta
tcell@meta.data$subcluster <- Idents(tcell)
tcell@meta.data$cell_name <- rownames(tcell@meta.data)
tcell_meta <- tcell@meta.data %>% select(cell_name, subcluster)
tcell_meta$tumor <- tcell_meta$subcluster

# B cell meta
Bcell@meta.data$subcluster <- Idents(Bcell)
Bcell@meta.data$cell_name <- rownames(Bcell@meta.data)
Bcell_meta <- Bcell@meta.data %>% select(cell_name, subcluster)
Bcell_meta$tumor <- Bcell_meta$subcluster

# Myeloid meta
Myeloi@meta.data$subcluster <- Idents(Myeloi)
Myeloi@meta.data$cell_name <- rownames(Myeloi@meta.data)
Myeloi_meta <- Myeloi@meta.data %>% select(cell_name, subcluster)
Myeloi_meta$tumor <- Myeloi_meta$subcluster

# Combine all metadata
meta <- rbind(nmf_meta, nomalignant_meta, stromal_meta, tcell_meta, Bcell_meta, Myeloi_meta)
allcel <- subset(allcell2, cells = rownames(meta))
saveRDS(allcel1, file.path("/hanlab/xueqinliu/test/00data/all_cell_ano.rds"))

# Merge metadata by cell_name
input_meta <- left_join(allcel@meta.data, meta, by = "cell_name")
allcel@meta.data$subcluster <- input_meta$subcluster
allcel@meta.data$tumor <- input_meta$tumor

# Select candidate cell types for cellphoneDB
candidate <- allcel@meta.data %>% 
  filter(tumor %in% c("malignant","FAP+ fibro","Tip-like ECs","Anti−inflammatory M",
                     "Pro-inflammatory M","neutrophil","Treg cells","Th 17 cells","HSP+ T cells"))
table(candidate$tumor)
allcel1 <- subset(allcel, cells = rownames(candidate))

#==================== Subset Intermediate group and timepoints ====================
Intermediate <- subset(allcel1, Group_Curativeeffect=="Intermediate")
Intermediate_CRC0 <- subset(Intermediate, TimePoints=="CRC_0")
Intermediate_CRC2 <- subset(Intermediate, TimePoints=="CRC_2")
Intermediate_CRC4 <- subset(Intermediate, TimePoints=="CRC_4")

#==================== Custom function: export h5ad + metadata txt for CellPhoneDB ====================
cpdb_export <- function(seu_obj, out_folder, suffix){
  setwd(out_folder)
  seu_fit <- UpdateSeuratObject(seu_obj)
  SaveH5Seurat(seu_fit, filename = paste0("Intermediate_fit_", suffix, ".h5Seurat"))
  Convert(paste0("Intermediate_fit_", suffix, ".h5Seurat"), dest = "h5ad")
  
  # prepare meta table for cellphoneDB
  seu_obj@meta.data$annotation <- seu_obj@meta.data$tumor
  meta_data <- data.frame(
    cell_names = rownames(seu_obj@meta.data),
    Cluster_type = seu_obj@meta.data$annotation
  )
  meta_data[is.na(meta_data)] = "Unknown" # Fix typo: Unkown -> Unknown
  write.table(meta_data, 
              paste0("Intermediate_CRC_", suffix, "_stromal_meta.txt"), 
              sep='\t', quote=F, row.names=F)
}

# Run export for CRC0, CRC2, CRC4
cpdb_export(Intermediate_CRC0, cpdb_crc0, "CRC0")
cpdb_export(Intermediate_CRC2, cpdb_crc2, "CRC2")
cpdb_export(Intermediate_CRC4, cpdb_crc4, "CRC4")