library(Seurat)
library(monocle)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

# ===================== Project Path Configuration =====================
# Root directory of the analysis project
proj_root <- "./test"
# Directory for raw seurat rds
raw_data_dir <- file.path(proj_root, "00data")
# Output directory for monocle results
monocle_out_dir <- file.path(proj_root, "02stromal", "05monocle")

# Create output folder automatically if it does not exist
if (!dir.exists(monocle_out_dir)) {
  dir.create(monocle_out_dir, recursive = TRUE)
}
# ======================================================================

# Load stromal Seurat object
stromal <- readRDS(file = file.path(raw_data_dir, "stromal_umap_new.rds"))
table(Idents(stromal))

# stromal <- subset(stromal, Group_Curativeeffect=="Resistant")
# stromal_Resistant <- subset(stromal, TimePoints=="CRC_Normal", invert=TRUE)

# Subset FAP+ fibro and remove CRC_Normal timepoint
FAP <- subset(stromal, idents = "FAP+ fibro")
FAP <- subset(FAP, TimePoints == "CRC_Normal", invert = TRUE)

# Extract count matrix for monocle CDS object
counts <- FAP@assays$RNA@counts
counts <- as.matrix(counts)

# Build phenoData: AnnotatedDataFrame, rows = cells, columns = cell metadata
pDat <- data.frame(cellname = colnames(FAP), FAP@meta.data)
pDat[is.na(pDat)] <- "Normal"

# Build featureData: AnnotatedDataFrame, rows = genes
fDat <- rownames(FAP)
fDat <- as.data.frame(fDat)
rownames(fDat) <- fDat$fDat
names(fDat) <- c("gene_short_name")

pd <- new("AnnotatedDataFrame", data = pDat)
fd <- new("AnnotatedDataFrame", data = fDat)

# Construct CDS object
HSMM3 <- newCellDataSet(as(counts, "sparseMatrix"),
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())

# Estimate size factors and dispersions
HSMM3 <- estimateSizeFactors(HSMM3)
HSMM3 <- estimateDispersions(HSMM3)
HSMM <- detectGenes(HSMM3, min_expr = 0.1)
print(head(fData(HSMM)))

## Quality control
# Filter expressed genes (detected in >=10 cells)
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
HSMM <- HSMM[expressed_genes, ]
print(head(pData(HSMM)))

# Filter cells with at least 200 expressed genes
valid_cells <- row.names(subset(pData(HSMM), num_genes_expressed >= 200))
length(valid_cells)
HSMM <- HSMM[, valid_cells]

expdt <- exprs(HSMM)
dim(expdt)

# Calculate total UMI counts per cell and add to metadata
pData(HSMM)$Total_mRNAs <- Matrix::colSums(expdt)
head(pData(HSMM))
dim(HSMM)

# Filter cells with total mRNA less than 1e6
valid_cells2 <- pData(HSMM)$Total_mRNAs < 1e6
HSMM <- HSMM[, valid_cells2]

# Filter cells within mean ± 2SD of log10(total mRNA)
upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))
valid_cells3 <- pData(HSMM)$Total_mRNAs > lower_bound & pData(HSMM)$Total_mRNAs < upper_bound
HSMM <- HSMM[, valid_cells3]
dim(HSMM)


## Trajectory construction
# Differential gene test for trajectory gene selection
clustering_DEG_genes <- differentialGeneTest(HSMM[expressed_genes,],
                                            fullModelFormulaStr = '~TimePoints',
                                            cores = 20)

head(clustering_DEG_genes)

# Select DEG with qval < 0.05 and sort by q-value
deg <- subset(clustering_DEG_genes, qval < 0.05)
deg <- deg[order(deg$qval, decreasing = F), ]
head(deg)

# Save DEG table
write.csv(deg, file = file.path(monocle_out_dir, "Resistant_subcluster_fap.DEG.csv"), row.names = F)
deg <- read.csv(file = file.path(monocle_out_dir, "Resistant_subcluster_fap.DEG.csv"), header = T)
rownames(deg) <- deg$gene_short_name
clustering_DEG_genes <- deg

# Select ordering genes
ordering_genes <- row.names(subset(clustering_DEG_genes, qval < 0.01))
HSMM_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:800]
filnaly <- Reduce(intersect, list(ordering_genes, HSMM_ordering_genes))

# Set ordering filter and visualize ordering genes
HSMM_myo <- setOrderingFilter(HSMM, ordering_genes = filnaly)
# pdf(file.path(monocle_out_dir, "01nores_monocle_CD4.ordergenes.pdf"))
plot_ordering_genes(HSMM_myo)
# dev.off()

# Dimension reduction and cell ordering by pseudotime
HSMM_myo <- reduceDimension(HSMM_myo, max_components = 2, method = 'DDRTree')
HSMM_myo <- orderCells(HSMM_myo)
# save(HSMM_myo, file = file.path(monocle_out_dir, "Resistant_subcluster_CRC0_HSMM_myo_malignant_normal.RData"))


cols_1 <- c("#008000","#0000FF","#DC143C","#BC8F8F")
cols_2 <- c("#f4a460","#ffd700","#90ee90","#48d1cc","#6495ed","#4169e1","#9370db","#d2691e","#db7093","#dcdcdc")

# Plot cell trajectory
pdf(file.path(monocle_out_dir, "01fap_psudotime.pdf"), width = 7, height = 7)
plot_cell_trajectory(HSMM_myo, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
plot_cell_trajectory(HSMM_myo, color_by = "State", size = 1, show_backbone = TRUE) + scale_colour_manual(values = cols_2)
plot_cell_trajectory(HSMM_myo, color_by = "TimePoints") + scale_colour_manual(values = cols_1)
plot_cell_trajectory(HSMM_myo, color_by = "TimePoints") + facet_wrap(~Group_Curativeeffect, nrow = 2) + scale_colour_manual(values = cols_1)
dev.off()


# Define color palette for violin plot
cols_1 <- c("#008000","#0000FF","#DC143C","#BC8F8F")

# Gene list definition
high_in_stage4 <- c("INHBA",  "HGF", "COL11A1")

# Custom function to generate violin plot with facet (monocle native plot cannot add facet_wrap)
make_gene_violin_facet <- function(gene, cds){
  meta_df <- pData(cds)
  meta_df$expr <- exprs(cds)[gene, ]
  ggplot(meta_df, aes(x = TimePoints, y = expr, fill = TimePoints)) +
    geom_violin(scale = "width") +
    geom_jitter(size = 0.3, alpha = 0.4) +
    facet_wrap(~Group_Curativeeffect, nrow = 3) +
    scale_fill_manual(values = cols_1) +
    theme_bw() +
    labs(title = gene, x = "TimePoints", y = "Expression")
}

# PDF output for violin plots
pdf(file.path(monocle_out_dir, "02fap_psudotime-gene.pdf"), width = 4, height = 4)
for(g in high_in_stage4){
  print(FeaturePlot(stromal, features = g))
  p <- make_gene_violin_facet(g, HSMM_myo)
  print(p)
}
dev.off()


# Test association between candidate genes and pseudotime using spline regression
candidate_genes <- read.csv(file = file.path("/hanlab/xueqinliu/test/02stromal", "01candidate_gene.csv")) ### shared genes in stage0&4 + stage4-specific upregulated genes compared to stage0 sensitive
candidate_genes <- candidate_genes$gene

Time_diff <- differentialGeneTest(HSMM_myo[candidate_genes,],
                                 cores = 1,
                                 fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[, c(5,2,3,4,1,6,7)] # Reorder columns
write.csv(Time_diff, file.path(monocle_out_dir, "Time_diff_all.csv"), row.names = F)
Time_genes <- Time_diff %>% pull(gene_short_name) %>% as.character()

p <- plot_pseudotime_heatmap(HSMM_myo[Time_genes,], num_clusters = 3, show_rownames = T, return_heatmap = T)
ggsave(file.path(monocle_out_dir, "02candidate_genes_heatmapAll.pdf"), p, width = 5, height = 15)


# BEAM test for branch-dependent gene expression at branch_point = 2
BEAM_res <- BEAM(HSMM_myo[candidate_genes], branch_point = 2, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval), ]
BEAM_res <- BEAM_res[, c("gene_short_name", "pval", "qval")]

pdf(file.path(monocle_out_dir, "03_3BEAMcandidate_genes_heatmapAll.pdf"), width = 5, height = 15)
plot_genes_branched_heatmap(HSMM_myo[row.names(subset(BEAM_res, qval < 1e-2)),],
                            branch_point = 2,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()