library(Seurat)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(tidyverse)
library(NMF)
library(DT)
library(ggsci)

# ===================== Path configuration: root directory, use file.path, remove all setwd =====================
root_dir <- "./test"
path_all_epi <- file.path(root_dir, "00data", "all_epi_nomalignant_malignant.rds")
out_nmf_main <- file.path(root_dir, "04epi", "01manlignant_NMF")
out_nmf_rank8 <- file.path(root_dir, "04epi", "01manlignant_NMF", "nmf_rank_8")
save_rds_path <- file.path(root_dir, "00data", "01_rank8_malignant.nmf.rds")

# Load seurat object and subset malignant cells
all_epi <- readRDS(path_all_epi)
malignant <- subset(all_epi, tomur == "malignant")
head(malignant@meta.data)

# Reconstruct seurat object, remove mitochondrial genes
malignant <- CreateSeuratObject(counts = malignant@assays$RNA@counts, meta.data = malignant@meta.data)
malignant[["percent.mt"]] <- PercentageFeatureSet(malignant, pattern = "^MT-")
mt.index <- grep(pattern = "^MT-", x = rownames(malignant@assays$RNA@counts), value = FALSE)
malignant <- malignant[-mt.index, ]

# Normalization
malignant <- NormalizeData(malignant, verbose = FALSE)
GetAssay(malignant, assay = "RNA")

# PCA workflow
malignant <- FindVariableFeatures(malignant, selection.method = 'vst')
var_genes <- VariableFeatures(malignant)
malignant <- ScaleData(malignant, features = var_genes, do.center = FALSE)
malignant <- RunPCA(object = malignant, pc.genes = var_genes, npcs = 80)

# Elbow plot for PCA
pdf(file.path(out_nmf_main, "ElbowPlot.pdf"), width = 6, height = 4)
print(ElbowPlot(malignant, ndims = 50))
dev.off()

# Clustering and UMAP based on PCA
malignant <- FindNeighbors(malignant, dims = 1:35)
res.used <- 0.23
malignant <- FindClusters(object = malignant, verbose = TRUE, resolution = res.used)
table(malignant$seurat_clusters)

set.seed(123)
malignant <- RunUMAP(malignant, reduction = "pca", dims = 1:35)

# Prepare scaled matrix as NMF input
vm <- malignant@assays$RNA@scale.data
dim(vm)

i0 <- which(colSums(vm) == 0)
i_na <- which(colSums(is.na(vm)) > 0)
vm_filtered <- vm[, -c(i0, i_na)]
dim(vm_filtered)

# Run NMF with rank=8
set.seed(123)
res <- nmf(vm, 8, nrun = 5, method = "snmf/r", seed = 'nndsvd')
save(res, file = file.path(out_nmf_main, "res8.RData"))

# WARNING: estim.coad object does not exist, original line commented
# pdf(file.path(out_nmf_main, "01_8nmf.pdf"))
# plot(estim.coad)
# dev.off()

# Load NMF result and assign NMF reduction to seurat object
load(file.path(out_nmf_main, "res8.RData"))
malignant@reductions$nmf <- malignant@reductions$pca
malignant@reductions$nmf@cell.embeddings <- t(coef(res))
malignant@reductions$nmf@feature.loadings <- basis(res)

# UMAP and clustering using NMF reduction
set.seed(219)
malignant.nmf <- RunUMAP(malignant, reduction = 'nmf', dims = 1:8) %>%
  FindNeighbors(reduction = 'nmf', dims = 1:8) %>%
  FindClusters()

# Assign cluster by max NMF coefficient
malignant.nmf$cluster <- apply(NMF::coefficients(res)[1:8, ], 2, which.max)

# Visualization of NMF clusters
pdf(file.path(out_nmf_rank8, "03malignant.nmf.pdf"))
print(DimPlot(malignant.nmf, group.by = "cluster") + scale_color_nejm())
print(DimPlot(malignant.nmf, group.by = "subcluster",
              cols = c("#f4a460","#ffd700","#90ee90","#48d1cc","#6495ed","#4169e1","#9370db","#d2691e","#db7093","#dcdcdc")))
dev.off()

# Extract top 300 marker features for each NMF component
f <- extractFeatures(res, 300L)
f <- lapply(f, function(x) rownames(res)[x])
f <- do.call("rbind", f)
DT::datatable(t(f))

# Export results to csv
write.csv(t(f), file.path(out_nmf_rank8, "top300.csv"), row.names = FALSE)
write.csv(table(malignant.nmf$orig.ident, malignant.nmf$cluster),
          file.path(out_nmf_rank8, "rank8_number.csv"), row.names = FALSE)

# DimPlot split by grouping variable
pdf(file.path(out_nmf_rank8, "04malignantnmf_groupby_subcluster.pdf"))
print(DimPlot(malignant.nmf,
              group.by = "subcluster",
              cols = c("#f4a460","#ffd700","#90ee90","#48d1cc","#6495ed","#4169e1","#9370db","#d2691e","#db7093","#dcdcdc"),
              split.by = "Group_Curativeeffect", ncol = 3) +
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")))
dev.off()

# Save final NMF seurat object
saveRDS(malignant.nmf, file = save_rds_path)