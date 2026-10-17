### Define malignant cell inferCNV heatmap
library(infercnv)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# ========= Set paths at the very beginning =========
setwd("./test/infercnv/result")

# Read infercnv object
infercnv_obj = readRDS("./test/infercnv/result/epi/C6/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data  # genes × cells

# Extract cell indices
normal_loc <- infercnv_obj@reference_grouped_cell_indices
normal_loc <- normal_loc$CRC_Normal
test_loc <- infercnv_obj@observation_grouped_cell_indices
test_loc0 <- test_loc$CRC_0
test_loc2 <- test_loc$CRC_2
test_loc4 <- test_loc$CRC_4

# Build cell annotation dataframe
anno.df=data.frame(
  CB = c(colnames(expr)[normal_loc],
         colnames(expr)[test_loc0],
         colnames(expr)[test_loc2],
         colnames(expr)[test_loc4]),
  class = c(rep("CRC_Normal",length(normal_loc)),
            rep("CRC_0",length(test_loc0)),
            rep("CRC_2",length(test_loc2)),
            rep("CRC_4",length(test_loc4)))
)

# Gene ordering file
gn <- rownames(expr)
geneFile <- read.table("/hanlab/xueqinliu/test/infercnv/input/Epithelial_C7_geneOrderingFile1.txt",
                      header = F,sep = "\t",stringsAsFactors = F)
# Chromosome vector with defined order
chr_order <- c(as.character(1:22),"X","Y")
geneFile <- filter(geneFile, V2 %in% chr_order)
rownames(geneFile)=geneFile$V1

# Keep intersect genes to synchronize genes and expression matrix
keep_genes <- intersect(gn, geneFile$V1)
sub_geneFile <-  geneFile[keep_genes,]
expr <- expr[keep_genes,]

# Build chromosome factor with fixed order
sub_geneFile$V2 <- factor(sub_geneFile$V2, levels = chr_order)

# ========= Kmeans clustering: cluster cells into 8 groups =========
set.seed(20210418)
# t(expr): cells × genes
kmeans.result <- kmeans(t(expr), centers = 8)
kmeans_df <- data.frame(kmeans_class = kmeans.result$cluster)
kmeans_df$CB = rownames(kmeans_df)
kmeans_df <- kmeans_df %>% inner_join(anno.df, by="CB") 

# Sort cells according to kmeans groups
kmeans_df_s <- arrange(kmeans_df, kmeans_class) 
rownames(kmeans_df_s) <- kmeans_df_s$CB
kmeans_df_s$CB <- NULL
kmeans_df_s$kmeans_class <- factor(kmeans_df_s$kmeans_class, levels = 1:8) 

# Color palette
color_v <- RColorBrewer::brewer.pal(11, "Spectral")[1:8] 
names(color_v) <- as.character(1:8)
anno_col <- list(
  class = c("CRC_Normal"="red","CRC_0" = "blue","CRC_2" = "green","CRC_4" = "orange"),
  kmeans_class = color_v
)

# ====== Top annotation: chromosomes (column annotation, genes grouped by chromosome) ======
# anno_block adds chromosome labels for column_split groups
top_anno <- HeatmapAnnotation(
  chr = anno_block(
    gp = gpar(fill = NA, col = NA), 
    labels = unique(sub_geneFile$V2), 
    labels_gp = gpar(cex = 0.8)
  ),
  show_legend = F
)

# Left row annotation (cell annotations)
left_anno <- rowAnnotation(df = kmeans_df_s, col = anno_col)

# Plot matrix: cell order matches annotation; t(expr) = cells as rows, genes as columns
plot_mat <- t(expr)[rownames(kmeans_df_s), ]
pdf("C6-1.pdf",width = 16,height = 10)
ht = Heatmap(
  plot_mat,
  # Color for inferCNV copy-number expression
  #col = colorRamp2(c(0.8,1,1.3), c("#377EB8","#F0F0F0","#E41A1C")), 
  col = colorRamp2(c(0.8,1,1.3), c("#377EB8","#F0F0F0","#E41A1C")), 
  cluster_rows = F,
  cluster_columns = F,
  show_column_names = F,
  show_row_names = F,
  
  # Split gene columns by chromosome
  column_split = sub_geneFile$V2, 
  column_gap = unit(2, "mm"),
  
  heatmap_legend_param = list(
    title = "Modified expression",
    direction = "vertical",
    title_position = "leftcenter-rot",
    at=c(0.8,1,1.3),
    legend_height = unit(3, "cm")
  ),
  
  top_annotation = top_anno,
  left_annotation = left_anno, 
  row_title = NULL,
  column_title = NULL
)
draw(ht, heatmap_legend_side = "right")
dev.off()