library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# ===================== Path configuration: root directory, use file.path for all paths, remove all setwd =====================
root_dir <- "./test"

# Input file paths
path_malignant_nmf <- file.path(root_dir, "00data", "01_rank8_malignant.nmf.rds")
path_all_epi <- file.path(root_dir, "00data", "all_epi_nomalignant_malignant.rds")
path_gene_csv <- file.path(root_dir, "04epi", "04manlignant_NMF", "nmf_rank_8", "81_up_gene.csv")

# Output pdf saved to corresponding subdirectory
out_dir_nmf <- file.path(root_dir, "04epi", "04manlignant_NMF", "nmf_rank_8")
out_dir_seurat_plot <- file.path(root_dir, "02zhoulian", "03seurat", "00_all")

# Read NMF malignant cell object
malignant <- readRDS(path_malignant_nmf)
Idents(malignant) <- malignant$subcluster

# UMAP plot 1
pdf(file.path(out_dir_seurat_plot, "01.pdf"), width = 4, height = 4)
DimPlot(malignant, reduction = "umap", pt.size = 0.5,
        cols = c("#74b4f3","#77d1fc","#f4cee1","#f6aba6","#9fd1a8","#f9deb5","#a3d1e1","#7cb87a","#fad7ce","#d592be"))
DimPlot(malignant, reduction = "umap", pt.size = 0.5,
        cols = c("#74b4f3","#77d1fc","#f4cee1","#f6aba6","#9fd1a8","#f9deb5","#a3d1e1","#7cb87a","#fad7ce","#d592be"),
        split.by = "Group_Curativeeffect")
dev.off()

# Read all_epi object
all_epi <- readRDS(path_all_epi)

# UMAP plot for all_epi
pdf(file.path(out_dir_seurat_plot, "02.pdf"), width = 4, height = 4)
DimPlot(all_epi, reduction = "umap", pt.size = 0.5, cols = c("#f6aba6","#77a1ec","#a3c4a5"))
dev.off()