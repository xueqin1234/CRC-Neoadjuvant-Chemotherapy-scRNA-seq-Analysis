# ===================== Load packages =====================
library(Seurat)
library(harmony)
library(cowplot)
library(ggplot2)
library(sctransform)
library(future)
library(glmGamPoi)
library(stringr)
library(DoubletFinder)

# ===================== Global parameters (modify here) =====================
# Root directory
root_dir <- "./test"
data_dir <- file.path(root_dir, "data")
out_dir <- file.path(root_dir, "01pre_deal_data")
ratio_dir <- file.path(out_dir, "cell_ratio")

# Rename cluster identities: named vector, old cluster name = new cell type
rename_idents_vec <- c(
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

# Cell types to remove
keep_cells <- c("Basal", "Acinar cells")
# Factor levels for cell annotation
cell_levels <- c("Epithelial cells","Stromal cells","T cells","B cells","Myeloid","Mast")
# Color palette for all cell types
cell_colors <- c("#FF3300","#9900FF","#6699FF","#009933","#00FF33","#FF9900")
# Levels and colors after removing epithelial cells
cell_levels_noEpi <- setdiff(cell_levels, "Epithelial cells")
cell_colors_noEpi <- cell_colors[cell_levels %in% cell_levels_noEpi]

# Marker genes for dot plot
markers.to.plot <- c("EPCAM","PHGR1","MUC13","ELF3","LGALS4","AGR2","ENG","PLVAP","TAGLN","RGS5","MMP2","COL1A2","COL3A1","THY1","CD3D","CD3E","CD3G","JCHAIN","BANK1","IGHG4","MZB1","S100A8","S100A9","LYZ","LST1","CD14","CCL3","IL1B","C1QB","KIT","TPSB2","CPA3","MS4A2")

# Sample order for plotting
sample_order <- c("C30","C40","C50","C60","C80","C100","C20","C70","C90","C110","C130","C120",
                  "C12","C32","C42","C52","C62","C102","C92","C122","C132","C44","C54","C64","C104",
                  "C94","C114","C134","N44","N64","N94","N114","N134")

# ===================== Function: generate stacked bar plot for cell proportion =====================
plot_cell_ratio <- function(seu_obj, out_csv, out_pdf, fill_colors, fill_levels, sample_order){
  # Count cells per cluster per sample
  tab <- table(Idents(seu_obj), seu_obj$orig.ident)
  # Normalize within each sample (column sum = 1)
  cell_proportion <- prop.table(tab, margin = 2)
  cell.prop <- as.data.frame(cell_proportion)
  colnames(cell.prop) <- c("cluster","sample","proportion")
  
  # Set factor levels to control plotting order
  cell.prop$cluster <- factor(cell.prop$cluster, levels = fill_levels)
  cell.prop$sample <- factor(cell.prop$sample, levels = sample_order)
  
  # Export proportion table
  write.csv(cell.prop, file.path(ratio_dir, out_csv), row.names = F)
  
  # ggplot for stacked bar
  p <- ggplot(cell.prop, aes(x = sample, y = proportion, fill = cluster)) +
    geom_col(position = "fill", width = 0.5) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(rep(1,4),'cm')) +
    scale_fill_manual(values = fill_colors) +
    scale_x_discrete(limits = sample_order) +
    coord_flip()
  
  ggsave(file.path(ratio_dir, out_pdf), plot = p, units="in", dpi=300, width=5, height=8, device="pdf")
  return(p)
}

# ===================== Main workflow =====================
# Create output directories if missing
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
if(!dir.exists(ratio_dir)) dir.create(ratio_dir, recursive = TRUE)

# Load Seurat object
allcell <- readRDS(file.path(data_dir, "03CRC_ALL_single_harmony.rds"))

# Rename original cluster ids to cell type names
allcell1 <- RenameIdents(allcell, rename_idents_vec)
table(Idents(allcell1))

# Subset: remove Basal and Acinar cells
allcell2 <- subset(allcell1, idents = keep_cells, invert = TRUE)
# Set cell type order
Idents(allcell2) <- factor(Idents(allcell2), levels = cell_levels)
allcell2$immune_annotation <- Idents(allcell2)
table(Idents(allcell2))

# -------------------- DotPlot for marker genes --------------------
pdf(file.path(out_dir, "02maker_CRC0.pdf"))
p1 <- DotPlot(allcell2, features = markers.to.plot, assay='RNA', group.by = 'immune_annotation') +
  theme(axis.text.x=element_text(face="italic",size=7, vjust = 1),
        axis.text.y = element_text(size=12)) +
  RotatedAxis() +
  coord_fixed(ratio=8/2)
print(p1)
dev.off()

# -------------------- UMAP visualization of cell types --------------------
pdf(file.path(out_dir, "01_CRC_cellnames.pdf"))
p_umap <- DimPlot(allcell2, reduction = "umap", pt.size = 0.5, cols = cell_colors)
print(p_umap)
dev.off()

# -------------------- Cell proportion plot (all retained cell types) --------------------
p_all <- plot_cell_ratio(seu_obj = allcell2,
                         out_csv = "sample_cell.prop.csv",
                         out_pdf = "02sample_cell_proportion.pdf",
                         fill_colors = cell_colors,
                         fill_levels = cell_levels,
                         sample_order = sample_order)

# -------------------- Cell proportion plot (remove epithelial cells) --------------------
allcell_noEpi <- subset(allcell2, idents = "Epithelial cells", invert = TRUE)
table(Idents(allcell_noEpi))

p_noEpi <- plot_cell_ratio(seu_obj = allcell_noEpi,
                           out_csv = "sample_cell_prop_noEpi.csv",
                           out_pdf = "02sample_to_evir_proportion.pdf",
                           fill_colors = cell_colors_noEpi,
                           fill_levels = cell_levels_noEpi,
                           sample_order = sample_order)