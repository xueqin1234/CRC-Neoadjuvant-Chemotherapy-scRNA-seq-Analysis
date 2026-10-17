#===================== Load Packages =====================
library(Seurat)
library(harmony)
library(cowplot)
library(ggplot2)
library(sctransform)
library(future)
library(glmGamPoi)
library(stringr)
library(DoubletFinder)
library(ggsci)
library(vcdExtra)
library(vcd)
library(pheatmap)

#===================== Project Path Config =====================
proj_root <- "./test"
raw_data_dir  <- file.path(proj_root, "data")
out_myeloid   <- file.path(proj_root, "05Myeloi")
out_roRE      <- file.path(out_myeloid, "01RORE")
out_data_dir  <- file.path(proj_root, "00data")

# Create output folders automatically
if (!dir.exists(out_myeloid)) dir.create(out_myeloid, recursive = TRUE)
if (!dir.exists(out_roRE))   dir.create(out_roRE, recursive = TRUE)

#===================== Load & Preprocess Seurat Object =====================
# Load full CRC object
allcell <- readRDS(file.path(raw_data_dir, "03CRC_ALL_single_harmony.rds"))

# Rename cluster identities
allcell1 <- RenameIdents(
  allcell,
  `0` = "T cells", `1` = "Epithelial cells", `2` = "T cells", `3` = "B cells",
  `4` = "Epithelial cells", `5` = "Myeloid", `6` = "Stromal cells", `7` = "T cells",
  `8` = "B cells", `9` = "B cells", `10` = "Epithelial cells", `11` = "Myeloid",
  `12` = "Stromal cells", `13` = "Epithelial cells", `14` = "Basal",`15` = "Stromal cells",
  `16` = "Mast",`17` = "Epithelial cells", `18` = "Myeloid",`19` = "Basal",
  `20` = "Epithelial cells",`21` = "B cells", `22` = "Basal",`23` = "Myeloid",
  `24` = "Epithelial cells",`25` = "Epithelial cells", `26` = "Epithelial cells",
  `27` = "T cells", `28` = "Stromal cells",`29` = "B cells", `30` = "Acinar cells",
  `31` = "Stromal cells", `32` = "Myeloid", `33` = "Epithelial cells", `34` = "T cells",
  `35` = "Epithelial cells", `36` = "T cells", `37` = "Stromal cells", `38` = "B cells",
  `39` = "B cells", `40` = "Myeloid"
)
table(Idents(allcell1))

# Filter: remove Basal, Acinar and C114 sample
allcell2 <- subset(allcell1, idents = c("Basal", "Acinar cells"), invert = TRUE)
allcell2 <- subset(allcell2, orig.ident == "C114", invert = TRUE)

# Set factor order for cell types
Idents(allcell2) <- factor(
  Idents(allcell2),
  levels = c("Epithelial cells","Stromal cells","T cells","B cells","Myeloid","Mast")
)
table(Idents(allcell2))

# Extract Myeloid + Mast
my  <- subset(allcell2, idents = "Myeloid")
ma  <- subset(allcell2, idents = "Mast")
cell_merge <- merge(my, ma)
Myeloi1 <- subset(allcell2, cells = rownames(cell_merge@meta.data))

#===================== Recreate Seurat Object & QC =====================
Myeloi1 <- CreateSeuratObject(counts = Myeloi1@assays$RNA@counts, meta.data = Myeloi1@meta.data)
Myeloi1[["percent.mt"]] <- PercentageFeatureSet(Myeloi1, pattern = "^MT-")
mt.index <- grep("^MT-", rownames(Myeloi1), value = FALSE)
Myeloi1 <- Myeloi1[-mt.index, ]

# Standard Seurat workflow
Myeloi1 <- NormalizeData(Myeloi1, verbose = FALSE)
Myeloi1 <- FindVariableFeatures(Myeloi1, selection.method = 'vst', nfeatures = 5000)
var_genes <- VariableFeatures(Myeloi1)
Myeloi1 <- ScaleData(Myeloi1, features = var_genes, verbose = FALSE)
Myeloi1 <- RunPCA(Myeloi1, features = var_genes, npcs = 80, verbose = FALSE)

# UMAP & Clustering
Myeloi1 <- RunUMAP(Myeloi1, reduction = "pca", dims = 1:40)
Myeloi1 <- FindNeighbors(Myeloi1, reduction = "pca", dims = 1:40)
Myeloi1 <- FindClusters(Myeloi1, resolution = 0.6)
Myeloi1

# Find all markers & save
markers_all <- FindAllMarkers(Myeloi1)
write.csv(markers_all, file.path(out_myeloid, "02_2000gene_dele_c114_noharmony_Myeloi1_cell_maker.csv"), row.names = F)

#===================== UMAP plots =====================
pdf(file.path(out_myeloid, "01_5000_06Myeloid_dele_c114_noharmony_TimePoints_Curativeeffect_orig.ident.pdf"))
DimPlot(Myeloi1, group.by = "TimePoints")
DimPlot(Myeloi1, group.by = "TWO_Curativeeffect")
DimPlot(Myeloi1, reduction = "umap", group.by = "orig.ident", pt.size = .1)
DimPlot(Myeloi1, reduction = "umap", label = TRUE, pt.size = .1)
dev.off()

#==================== Add metadata: Group_Curativeeffect & Curativeeffect ====================
# Helper function to extract sample ID from cell name
get_sample_id <- function(cell_name){
  unlist(strsplit(cell_name, "\\_"))[1]
}
sample_vec <- sapply(colnames(Myeloi1), get_sample_id)

# Mapping table, easier to maintain than nested ifelse
map_group_cure <- c(
  C100="Resistant",C102="Resistant",C104="Resistant",C110="Sensitive",C12="Resistant",
  C120="Unknown",C122="Sensitive",C130="Sensitive",C132="Sensitive",C134="Sensitive",
  C20="Sensitive",C30="Intermediate",C32="Intermediate",C40="Intermediate",C42="Intermediate",
  C44="Intermediate",C50="Resistant",C52="Resistant",C54="Resistant",C60="Resistant",
  C62="Resistant",C64="Resistant",C70="Sensitive",C80="Resistant",C90="Sensitive",
  C92="Sensitive",C94="Sensitive",N114="Sensitive",N134="Sensitive",N44="Intermediate",N64="Resistant"
)
# unmatched default Sensitive
Group_Curativeeffect <- ifelse(sample_vec %in% names(map_group_cure), map_group_cure[sample_vec], "Sensitive")
Myeloi1 <- AddMetaData(Myeloi1, metadata = Group_Curativeeffect, col.name = "Group_Curativeeffect")

map_cure <- c(
  C20="Sensitive",C30="Sensitive",C40="Sensitive",C50="Resistant",C60="Resistant",
  C70="Sensitive",C90="Sensitive",C100="Resistant",C110="Sensitive",C120="Unknown",
  C130="Sensitive",C32="Sensitive",C12="Sensitive",C42="Sensitive",C62="Resistant",
  C92="Sensitive",C102="Resistant",C112="Sensitive",C132="Sensitive",C44="Sensitive",
  C54="Resistant",C64="Resistant",C94="Sensitive",C104="Resistant",C34="Sensitive",
  N44="Sensitive",C52="Resistant",C80="Resistant",C114="Sensitive",N114="Sensitive",
  C134="Sensitive",N134="Sensitive",N64="Resistant"
)
Curativeeffect <- ifelse(sample_vec %in% names(map_cure), map_cure[sample_vec], "Sensitive")
Myeloi1 <- AddMetaData(Myeloi1, metadata = Curativeeffect, col.name = "Curativeeffect")

Myeloi2 <- Myeloi1

#==================== Rename clusters & filter cell types ====================
Myeloi2 <- RenameIdents(
  Myeloi2,
  `0` = "neutrophil", `1` = "neutrophil", `2` = "mono-like", `3` = "Mast",
  `4` = "neutrophil", `5` = "Pro-inflammatory M", `6` = "epi", `7` = "Anti−inflammatory M",
  `8` = "neutrophil", `9` = "Mast", `10` = "neutrophil", `11` = "T", `12` = "Conventional DC2",
  `13` = "neutrophil", `14` = "epi",`15` = "LYVE1+ FOLR2+ macrophages", `16` = "B",
  `17` = "Conventional DC1", `18` = "Conventional DC1",`19` = "Plasmacytoid DC",
  `20` = "B", `21` = "Conventional DC1",`22` = "epi",`23` = "Mast",`24` = "Mast"
)
table(Idents(Myeloi2))

# Remove epi, T, B; retain target myeloid populations
keep_celltypes <- c("mono-like","Anti−inflammatory M","LYVE1+ FOLR2+ macrophages",
                    "Pro-inflammatory M","Conventional DC1","Conventional DC2",
                    "Plasmacytoid DC","neutrophil","Mast")
Myeloi2 <- subset(Myeloi2, idents = keep_celltypes)
# Remove sample C120
Myeloi2 <- subset(Myeloi2, orig.ident == "C120", invert = TRUE)

# Set factor order
Idents(Myeloi2) <- factor(
  Idents(Myeloi2),
  levels = keep_celltypes
)
table(Idents(Myeloi2))
Myeloi2@meta.data$cluster_annotation <- Idents(Myeloi2)

#==================== UMAP with annotated cell type ====================
color_myeloid <- c("#ffcc33","#cc9933","#996633","#9900ff","#ffa07a","#9999cc","#99ffff","#6600ff","#33ff33","#339933","#cc3300")
pdf(file.path(out_myeloid, "01new_myeloid_cellnames.pdf"))
DimPlot(Myeloi2, reduction = "umap", pt.size = 0.5, cols = color_myeloid)
dev.off()

#==================== Dotplot marker genes ====================
markers.to.plot <- c("CD14","FCGR3A","S100A10","LYZ","CD68","APOC1","APOE","SELENOP",
                     "C1QA","C1QB","SPP1","CD163","MRC1","LYVE1","FOLR2","RNASE1",
                     "IL1B","EREG","TIMP1","CD86","IDO1","BATF3","BIRC3","CLEC10A",
                     "CD1C","FCER1A","IRF8","CLEC4C","IL3RA","FCGR3B","IFITM2",
                     "CXCR2","CSF3R","S100A8","S100A9","TPSAB1","CPA3","MS4A2")

pdf(file.path(out_myeloid, "02new—maker_Myeloi2.pdf"), width = 12, height = 5)
DotPlot(Myeloi2, features = markers.to.plot, assay = 'RNA', group.by = 'cluster_annotation') +
  RotatedAxis() +
  theme(axis.text = element_text(size=10),
        axis.text.x = element_text(vjust = 1, size=12),
        axis.text.y = element_text(size=14))
dev.off()

#==================== Cell proportion per sample ====================
cell.tab <- table(Idents(Myeloi2), Myeloi2$orig.ident)
cell.prop.df <- as.data.frame(prop.table(cell.tab, margin = 2))
colnames(cell.prop.df) <- c("cluster","sample","proportion")

write.csv(cell.prop.df, file.path(out_myeloid, "1sample_cell.prop.csv"), row.names = F)
write.csv(as.data.frame(cell.tab), file.path(out_myeloid, "1sample_cell.prop2.csv"), row.names = F)

sample_order <- c("C20","C70","C90","C110","C130","C30","C40","C80","C50","C60","C100",
                  "C92","C122","C132","C32","C42","C52","C62","C102","C12","C94","C134",
                  "C44","C54","C64","C104","N44","N64","N94","N114","N134")

p_bar <- ggplot(cell.prop.df, aes(x = sample, y = proportion, fill = cluster)) +
  geom_col(position = "fill", width = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = color_myeloid) +
  scale_x_discrete(limits = sample_order) +
  theme(plot.margin = unit(rep(1,4), 'cm'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(file.path(out_myeloid, "03Myeloi_per_sample_cell_proportion.pdf"),
       plot = p_bar, units="in", dpi=300, width=10, height=3, device="pdf")

#==================== Contingency table & heatmap (RORE) ====================
# Subset by TimePoints, remove Unknown
Myeloi_CRC_0 <- subset(Myeloi2, TimePoints=="CRC_0")
Myeloi_CRC_0 <- subset(Myeloi_CRC_0, Group_Curativeeffect!="Unknown")
tab0 <- table(Myeloi_CRC_0$cluster_annotation, Myeloi_CRC_0$Group_Curativeeffect)
colnames(tab0) <- c("Intermediate_0","Resistant_0","Sensitive_0")

Myeloi_CRC_2 <- subset(Myeloi2, TimePoints=="CRC_2")
Myeloi_CRC_2 <- subset(Myeloi_CRC_2, Group_Curativeeffect!="Unknown")
tab2 <- table(Myeloi_CRC_2$cluster_annotation, Myeloi_CRC_2$Group_Curativeeffect)
colnames(tab2) <- c("Intermediate_2","Resistant_2","Sensitive_2")

Myeloi_CRC_4 <- subset(Myeloi2, TimePoints=="CRC_4")
Myeloi_CRC_4 <- subset(Myeloi_CRC_4, Group_Curativeeffect!="Unknown")
tab4 <- table(Myeloi_CRC_4$cluster_annotation, Myeloi_CRC_4$Group_Curativeeffect)
colnames(tab4) <- c("Intermediate_4","Resistant_4","Sensitive_4")

Myeloi_CRC_Normal <- subset(Myeloi2, TimePoints=="CRC_Normal")
tabN <- table(Myeloi_CRC_Normal$cluster_annotation, Myeloi_CRC_Normal$TimePoints)
colnames(tabN) <- "CRC_Normal"

specificity_mat <- cbind(tab0, tab2, tab4, tabN)

# CMH test & obs/exp ratio
CMHtest(specificity_mat)
chi_res <- chisq.test(specificity_mat)
expected_mat <- chi_res$expected
ratio_mat <- specificity_mat / expected_mat

df_scaled <- t(scale(t(ratio_mat), center = TRUE, scale = TRUE))
write.csv(df_scaled, file.path(out_roRE, "Tcell_all_R_obs_exp.csv"))

# Annotation for pheatmap
ann_col <- data.frame(
  TimePoints = factor(rep(c("CRC_0","CRC_2","CRC_4","Normal"), c(3,3,3,1))),
  Group = factor(c(rep(c("Intermediate","Resistant","Sensitive"),3), "Normal"))
)
rownames(ann_col) <- colnames(df_scaled)

ann_colors <- list(
  TimePoints = c(CRC_0 = "#008000", CRC_2 = "#0000FF", CRC_4 = "#DC143C", Normal="#BC8F8F"),
  Group = c(Intermediate="#2F4F4F",Resistant="#B8860B",Sensitive="#008B8B",Normal="#BC8F8F")
)

bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pdf(file.path(out_roRE, "01_new01Tcell_Group_Curativeeffect_R_obs_exp.pdf"), width=5, height=4)
pheatmap(df_scaled,
         color = c(colorRampPalette(c("#0000cc","#f5f5f5"))(length(bk)/2),
                   colorRampPalette(c("#f5f5f5","#cc0000"))(length(bk)/2)),
         breaks=bk,legend_breaks=seq(-2,2,2),
         cluster_row = FALSE, cluster_col = FALSE,
         annotation_col = ann_col, annotation_colors = ann_colors)
pheatmap(df_scaled, color = colorRampPalette(c("#0000cc","#f5f5f5","#cc6600","#cc0000"))(50),
         cluster_row = FALSE, cluster_col = FALSE, annotation_col = ann_col, annotation_colors = ann_colors)
pheatmap(df_scaled, color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_row = FALSE, cluster_col = FALSE, annotation_col = ann_col, annotation_colors = ann_colors)
dev.off()

# Save final myeloid object
saveRDS(Myeloi2, file.path(out_data_dir, "01Myeloi_umap.rds"))