library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)

# Define the root directory. All subsequent paths are built from this, no more setwd
root_dir <- "./test"

# Read CNV results
merge.CNV <- read.csv(
  file.path(root_dir, "infercnv", "result.csv"),
  check.names = FALSE
)
tumor.cell.names <- colnames(merge.CNV)
length(tumor.cell.names)

# Read epithelial Seurat object
all_epi <- readRDS(file.path(root_dir, "00data", "all_epi_umap.rds"))

head(all_epi@meta.data)

# Subset malignant cells, exclude CRC_Normal
malignant <- subset(all_epi, cells = tumor.cell.names)
malignant <- subset(malignant, TimePoints == "CRC_Normal", invert = TRUE)
malignant_meta <- malignant@meta.data
malignant_meta$tomur <- "malignant"

# Normal epithelial cells
normal_epi <- subset(all_epi, TimePoints == "CRC_Normal")
normal_epi_meta <- normal_epi@meta.data
normal_epi_meta$tomur <- "normal_epi"

# Non-malignant epithelial cells
nomalignant <- subset(all_epi, cells = tumor.cell.names, invert = TRUE)
nomalignant <- subset(nomalignant, TimePoints == "CRC_Normal", invert = TRUE)
nomalignant_meta <- nomalignant@meta.data
nomalignant_meta$tomur <- "nomalignant"

# Merge metadata, bind_rows is more stable
meta <- bind_rows(malignant_meta, normal_epi_meta, nomalignant_meta)
meta <- meta[, c("cell_name", "tomur")]

all_epi_meta <- all_epi@meta.data
all_epi_meta$tomur <- NULL

input_meta <- left_join(all_epi_meta, meta, by = "cell_name")
identical(input_meta$cell_name, rownames(all_epi@meta.data))
all_epi@meta.data$tomur <- input_meta$tomur

# Save output object
save_path <- file.path(root_dir, "00data", "all_epi_nomalignant_malignant.rds")
saveRDS(all_epi, file = save_path)

# Optional: check cell counts for each tomur group
table(all_epi$tomur)