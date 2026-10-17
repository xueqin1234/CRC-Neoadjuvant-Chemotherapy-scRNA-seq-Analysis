# ===================== Load packages =====================
library(Seurat)
library(ggplot2)
library(ggsci)
library(vcdExtra)
library(vcd)
library(pheatmap)

# ===================== Global parameters (modify here) =====================
root_dir <- "./test"
data_dir <- file.path(root_dir, "data")
stromal_out_dir <- file.path(root_dir, "02stromal")
rore_dir <- file.path(stromal_out_dir, "01rore")
save_rds_path <- file.path(root_dir, "00data/stromal_umap_new.rds")

# Cell type colors for stromal subclusters
stromal_colors <- c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF",
                    "#197EC0FF","#46732EFF","#71D0F5FF","#F79D1EFF","#1A9993FF","#FD8CC1FF","#7876B1FF")

# Sample order for proportion bar plot
sample_order <- c("C120","C20","C90","C92","C94","C110","C122","C130","C132","C134",
                  "C30","C32","C40","C42","C44","C70","C80","C50","C52","C54",
                  "C60","C62","C64","C100","C102","C104","C12","N44","N64","N94","N114","N134")

# Stromal subcluster rename vector
rename_stroma_vec <- c(
  `0` = "FAP+ fibro", `1` = "CCL8 fibro", `2` = "Unknown", `3` = "Pericyte",
  `4` = "Capillary ECs", `5` = "Tip-like ECs", `6` = "CCL8 fibro", `7` = "FAP+ fibro",
  `8` = "Normal mucosa fibro", `9` = "Venous ECs", `10` = "BMP fibro", `11` = "Myofibroblasts",
  `12` = "T", `13` = "Pericyte", `14` = "Enteric glial cells",`15` = "Smooth muscle cells",
  `16` = "Normal mucosa fibro",`17` = "BMP fibro", `18` = "CCL8 fibro",`19` = "B",
  `20` = "FAP+ fibro",`21` = "Lymphatic ECs", `22` = "Smooth muscle cells",`23` = "epi"
)

# Keep stromal subclusters after renaming
keep_stroma_subtypes <- c("FAP+ fibro","Myofibroblasts","CCL8 fibro","BMP fibro",
                          "Normal mucosa fibro","Tip-like ECs","Capillary ECs","Venous ECs",
                          "Lymphatic ECs","Pericyte","Enteric glial cells","Smooth muscle cells","Unknown")

# Factor level order for stromal subtypes
stroma_levels <- keep_stroma_subtypes

# Annotation config for pheatmap
anno_col_time <- factor(rep(c("CRC_0","CRC_2","CRC_4","Normal"), c(3, 3, 3,1)))
anno_col_group <- factor(c(rep(c("Intermediate","Resistant","Sensitive"),3), "Normal"))
ann_colors <- list(
  TimePoints = c(CRC_0 = "#008000", CRC_2 = "#0000FF", CRC_4 = "#DC143C","Normal"="#BC8F8F"),
  Group = c(Intermediate="#2F4F4F",Resistant="#B8860B",Sensitive="#008B8B",Normal="#BC8F8F")
)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))

# ===================== Function: plot cell proportion stacked bar =====================
plot_stroma_prop <- function(seu_obj, out_csv, out_pdf, fill_colors, sample_order){
  tab <- table(Idents(seu_obj), seu_obj$orig.ident)
  cell_proportion <- prop.table(tab, margin = 2)
  cell.prop <- as.data.frame(cell_proportion)
  colnames(cell.prop) <- c("cluster","sample","proportion")
  
  cell.prop$cluster <- factor(cell.prop$cluster, levels = stroma_levels)
  cell.prop$sample <- factor(cell.prop$sample, levels = sample_order)
  
  write.csv(cell.prop, file.path(stromal_out_dir, out_csv), row.names = F)
  write.csv(as.data.frame(cell_proportion), file.path(stromal_out_dir, "1sample_cell.prop2.csv"), row.names = T)
  
  p <- ggplot(cell.prop, aes(x = sample, y = proportion, fill = cluster)) +
    geom_col(position = "fill", width = 0.5) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(rep(1,4),'cm')) +
    scale_fill_manual(values = fill_colors) +
    scale_x_discrete(limits = sample_order) +
    coord_flip()
  
  ggsave(file.path(stromal_out_dir, out_pdf), plot = p, units="in", dpi=300, width=5, height=8, device="pdf")
  return(p)
}

# ===================== Main workflow =====================
# Create folders
if(!dir.exists(stromal_out_dir)) dir.create(stromal_out_dir, recursive = TRUE)
if(!dir.exists(rore_dir)) dir.create(rore_dir, recursive = TRUE)
if(!dir.exists(dirname(save_rds_path))) dir.create(dirname(save_rds_path), recursive = TRUE)

# Load full seurat object
setwd(data_dir)
allcell <- readRDS("03CRC_ALL_single_harmony.rds")

# Rename original clusters
allcell1 <- RenameIdents(allcell,
  `0` = "T cells", `1` = "Epithelial cells", `2` = "T cells", `3` = "B cells",
  `4` = "Epithelial cells", `5` = "Myeloid", `6` = "Stromal cells", `7` = "T cells",
  `8` = "B cells", `9` = "B cells", `10` = "Epithelial cells", `11` = "Myeloid",
  `12` = "Stromal cells", `13` = "Epithelial cells", `14` = "Basal",`15` = "Stromal cells",
  `16` = "Mast",`17` = "Epithelial cells", `18` = "Myeloid",`19` = "Basal", `20` = "Epithelial cells",
  `21` = "B cells", `22` = "Basal",`23` = "Myeloid", `24` = "Epithelial cells",`25` = "Epithelial cells",
  `26` = "Epithelial cells",`27` = "T cells", `28` = "Stromal cells",`29` = "B cells", `30` = "Acinar cells",
  `31` = "Stromal cells", `32` = "Myeloid", `33` = "Epithelial cells", `34` = "T cells",
  `35` = "Epithelial cells", `36` = "T cells", `37` = "Stromal cells", `38` = "B cells",
  `39` = "B cells", `40` = "Myeloid"
)
table(Idents(allcell1))

# Filter: remove Basal, Acinar cells and C114 sample
allcell2 <- subset(allcell1, idents = c("Basal", "Acinar cells"), invert = TRUE)
allcell2 <- subset(allcell2, orig.ident == "C114", invert = TRUE)

Idents(allcell2) <- factor(Idents(allcell2), levels = c("Epithelial cells","Stromal cells","T cells","B cells","Myeloid","Mast"))
table(Idents(allcell2))

# Subset stromal cells only
Stroma1 <- subset(allcell2, idents = "Stromal cells")

# Re-create new Seurat object from stromal raw count
Stroma <- CreateSeuratObject(counts = Stroma1@assays$RNA@counts, meta.data = Stroma1@meta.data)

# Remove mitochondrial genes
mt.index <- grep("^MT-", rownames(Stroma), value = FALSE)
if(length(mt.index) > 0){
  Stroma <- Stroma[-mt.index, ]
}

# Standard preprocessing
Stroma <- NormalizeData(Stroma, verbose = FALSE)
Stroma <- FindVariableFeatures(Stroma, selection.method='vst', nfeatures = 5000)
var_genes <- VariableFeatures(Stroma)
Stroma <- ScaleData(Stroma, features = var_genes, verbose = FALSE)
Stroma <- RunPCA(Stroma, features = var_genes, npcs = 80, verbose = FALSE)

# Elbow plot
pdf(file.path(stromal_out_dir, "elbow_stroma.pdf"))
ElbowPlot(Stroma, ndims = 50)
dev.off()

# UMAP & clustering
Stroma <- RunUMAP(Stroma, reduction = "pca", dims = 1:40, verbose = FALSE)
Stroma <- FindNeighbors(Stroma, reduction = "pca", dims = 1:40, verbose = FALSE)
Stroma <- FindClusters(Stroma, resolution = 0.5, verbose = FALSE)

# ---------------- Metadata: add treatment response groups ----------------
# Extract sample ID from cell name
get_sample_id <- function(cell_name){
  unlist(strsplit(cell_name, "\\_"))[1]
}
sample_vec <- sapply(colnames(Stroma), get_sample_id)
table(sample_vec)

# Group_Curativeeffect mapping
map_group1 <- c(
  C100="Resistant",C102="Resistant",C104="Resistant",C110="Sensitive",C12="Resistant",
  C120="Unknown",C122="Sensitive",C130="Sensitive",C132="Sensitive",C134="Sensitive",
  C20="Sensitive",C30="Intermediate",C32="Intermediate",C40="Intermediate",C42="Intermediate",
  C44="Intermediate",C50="Resistant",C52="Resistant",C54="Resistant",C60="Resistant",
  C62="Resistant",C64="Resistant",C70="Sensitive",C80="Intermediate",C90="Sensitive",
  C92="Sensitive",C94="Sensitive",N114="Sensitive",N134="Sensitive",N44="Intermediate",N64="Resistant"
)
Stroma$Group_Curativeeffect <- unname(map_group1[sample_vec])

# Curativeeffect mapping
map_group2 <- c(
  C20="Sensitive",C30="Sensitive",C40="Sensitive",C50="Resistant",C60="Resistant",
  C70="Sensitive",C90="Sensitive",C100="Resistant",C110="Sensitive",C120="Unknown",
  C130="Sensitive",C32="Sensitive",C12="Sensitive",C42="Sensitive",C62="Resistant",
  C92="Sensitive",C102="Resistant",C112="Sensitive",C132="Sensitive",C44="Sensitive",
  C54="Resistant",C64="Resistant",C94="Sensitive",C104="Resistant",C34="Sensitive",
  N44="Sensitive",C52="Resistant",C80="Resistant",C114="Sensitive",N114="Sensitive",
  C134="Sensitive",N134="Sensitive",N64="Resistant"
)
Stroma$Curativeeffect <- unname(map_group2[sample_vec])
table(Stroma$Curativeeffect)

# Rename stromal clusters to cell subtypes
Stroma1 <- RenameIdents(Stroma, rename_stroma_vec)

# Remove unwanted cell types (T, B, epi)
Stroma1 <- subset(Stroma1, idents = keep_stroma_subtypes)
Idents(Stroma1) <- factor(Idents(Stroma1), levels = stroma_levels)
table(Idents(Stroma1))
Stroma1$subcluster <- Idents(Stroma1)

# ---------------- UMAP plot of stromal subtypes ----------------
pdf(file.path(stromal_out_dir,"04sub_Stromal_cellnames.pdf"), width=8, height=5)
DimPlot(Stroma1, reduction = "umap", pt.size = 0.5, cols = stromal_colors)
dev.off()

# ---------------- Contingency table & CMH test + Obs/Exp heatmap ----------------
# Subset each timepoint and exclude Unknown group
stromal_CRC_0 <- subset(Stroma1, TimePoints=="CRC_0")
stromal_CRC_0 <- subset(stromal_CRC_0, Group_Curativeeffect!="Unknown")
spec_CRC0 <- table(stromal_CRC_0$subcluster, stromal_CRC_0$Group_Curativeeffect)
colnames(spec_CRC0) <- c("Intermediate_0","Resistant_0","Sensitive_0")

stromal_CRC_2 <- subset(Stroma1, TimePoints=="CRC_2")
stromal_CRC_2 <- subset(stromal_CRC_2, Group_Curativeeffect!="Unknown")
spec_CRC2 <- table(stromal_CRC_2$subcluster, stromal_CRC_2$Group_Curativeeffect)
colnames(spec_CRC2) <- c("Intermediate_2","Resistant_2","Sensitive_2")

stromal_CRC_4 <- subset(Stroma1, TimePoints=="CRC_4")
stromal_CRC_4 <- subset(stromal_CRC_4, Group_Curativeeffect!="Unknown")
spec_CRC4 <- table(stromal_CRC_4$subcluster, stromal_CRC_4$Group_Curativeeffect)
colnames(spec_CRC4) <- c("Intermediate_4","Resistant_4","Sensitive_4")

stromal_CRC_Normal <- subset(Stroma1, TimePoints=="CRC_Normal")
spec_Normal <- table(stromal_CRC_Normal$subcluster, stromal_CRC_Normal$TimePoints)

# Combine matrix
spec_all <- cbind(spec_CRC0, spec_CRC2, spec_CRC4, spec_Normal)
CMH_result <- CMHtest(spec_all)
print(CMH_result)

# Observed / Expected ratio and scale
chi_res <- chisq.test(spec_all)
exp_mat <- chi_res$expected
obs_exp_ratio <- spec_all / exp_mat

df_t <- t(obs_exp_ratio)
scaled_mat <- scale(df_t, center = TRUE, scale = TRUE)
all_mat <- t(scaled_mat)

# Export matrix
write.csv(all_mat, file.path(rore_dir,"stromal_all_R_obs_exp.csv"), row.names = TRUE)

# Prepare annotation for heatmap
annotation_col <- data.frame(TimePoints = anno_col_time, Group = anno_col_group)
rownames(annotation_col) <- colnames(all_mat)

# Plot pheatmap
pdf(file.path(rore_dir,"01_new_stromal_Group_Curativeeffect_R_obs_exp.pdf"), width=5, height=4)
pheatmap(all_mat,
         color = c(colorRampPalette(c("#0000cc","#f5f5f5"))(length(bk)/2),
                   colorRampPalette(c("#f5f5f5","#cc0000"))(length(bk)/2)),
         breaks=bk, legend_breaks=seq(-2,2,2),
         cluster_row = FALSE, cluster_col = FALSE,
         annotation_col = annotation_col, annotation_colors = ann_colors)
pheatmap(all_mat, color = colorRampPalette(c("#0000cc","#f5f5f5","#cc6600","#cc0000"))(50),
         cluster_row = FALSE, cluster_col = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors)
pheatmap(all_mat, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_row = FALSE, cluster_col = FALSE, annotation_col = annotation_col, annotation_colors = ann_colors)
dev.off()

# ---------------- Sample-wise stromal subcluster proportion barplot ----------------
p_stroma_prop <- plot_stroma_prop(seu_obj = Stroma1,
                                  out_csv = "1sample_cell.prop.csv",
                                  out_pdf = "03stromal_sample_cell_proportion.pdf",
                                  fill_colors = stromal_colors,
                                  sample_order = sample_order)

# Save final stromal seurat object
saveRDS(Stroma1, file = save_rds_path)

