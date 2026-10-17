#==================== Load Packages ====================
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
library(dplyr)

#==================== 【Path Configuration, only modify here】 ====================
root_dir <- "/hanlab/xueqinliu/test/04Bcell"
# Subdirectories
dir_04Bcell  <- file.path(root_dir, "04Bcell")
dir_ratio    <- file.path(dir_04Bcell, "01ratio")
dir_rore     <- file.path(dir_04Bcell, "RORE")

# Automatically create folders
if (!dir.exists(root_dir)) dir.create(root_dir, recursive = TRUE)
if (!dir.exists(dir_04Bcell)) dir.create(dir_04Bcell, recursive = TRUE)
if (!dir.exists(dir_ratio)) dir.create(dir_ratio, recursive = TRUE)
if (!dir.exists(dir_rore)) dir.create(dir_rore, recursive = TRUE)

# Path to raw RDS file (read original data, do not modify)
data_raw_path <- "/hanlab/xueqin/03zhoulian/data/03CRC_ALL_single_harmony.rds"

#==================== 1. Load raw object, extract B cell subset and re-analysis ====================
allcell <- readRDS(data_raw_path)

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
# Remove Basal and Acinar cells; remove sample C114
allcell2 <- subset(allcell1, idents = c("Basal", "Acinar cells"), invert = TRUE)
allcell2 <- subset(allcell2, orig.ident == "C114", invert = TRUE)

Idents(allcell2) <- factor(Idents(allcell2), 
                           levels = c("Epithelial cells","Stromal cells","T cells","B cells","Myeloid","Mast"))
table(Idents(allcell2))

# Subset B cells
Bcell <- subset(allcell2, idents = "B cells")

# Create new Seurat object, remove mitochondrial genes
Bcell <- CreateSeuratObject(counts = Bcell@assays$RNA@counts, meta.data = Bcell@meta.data)
Bcell[["percent.mt"]] <- PercentageFeatureSet(Bcell, pattern = "^MT-")
mt.index <- grep("^MT-", rownames(Bcell), value = FALSE)
Bcell <- Bcell[-mt.index, ]

Bcell <- NormalizeData(Bcell, verbose=FALSE)
Bcell <- FindVariableFeatures(Bcell, selection.method='vst', nfeatures = 2000)
gene <- VariableFeatures(Bcell)
Bcell <- ScaleData(Bcell, features=gene, verbose=FALSE) 
Bcell <- RunPCA(Bcell, npcs = 80, verbose=FALSE)

# Optional:
# DimHeatmap(Bcell, dims = 1:50, balanced = TRUE)
# ElbowPlot(Bcell, ndims = 50)

Bcell <- RunUMAP(Bcell, reduction = "pca", dims = 1:40, verbose=FALSE)
Bcell <- FindNeighbors(Bcell, reduction = "pca", dims = 1:40, verbose=FALSE)
Bcell <- FindClusters(Bcell, resolution = 0.6, verbose=FALSE)

# Find marker genes and export CSV to 04Bcell folder
setwd(dir_04Bcell)
df <- FindAllMarkers(Bcell)
write.csv(df, file="02_2000_0.6gene_dele_c114_noharmony_Bcell_cell_maker.csv", row.names = F)

# UMAP plot of raw clustering
pdf("0102_2000_0.6_Bcell_cell_maker.pdf")
DimPlot(Bcell, reduction = "umap", label = TRUE, pt.size = .1)
dev.off() 

# Rename B cell clusters
Bcell1 <- RenameIdents(Bcell, 
                       `0` = "MALT B cells", 
                       `1` = "Follicular B cells", 
                       `2` = "MALT B cells", 
                       `3` = "MALT B cells", 
                       `4` = "Naive B cells", 
                       `5` = "Unknown", 
                       `6` ="Tumor-associated B cells", 
                       `7` = "Plasma IgG", 
                       `8` = "Follicular B cells", 
                       `9` = "MALT B cells", 
                       `10` = "MALT B cells", 
                       `11` = "MALT B cells", 
                       `12` = "Tumor-associated B cells", 
                       `13` = "Plasma IgG", 
                       `14` = "T",
                       `15` = "Plasma IgG prolif", 
                       `16` = "Germinal center B cells",
                       `17` = "MALT B cells", 
                       `18` = "Plasma IgG",
                       `19` = "MALT B cells",
                       `20` = "Naive B cells", 
                       `21` = "MALT B cells",
                       `22` = "Tumor-associated B cells",
                       `23` = "Plasma IgG",
                       `24` = "inflammation-associated B cells")
table(Idents(Bcell1))

# Remove contaminated T cells, retain target B cell subsets
Bcell1 <- subset(Bcell1, idents = c("Naive B cells","Follicular B cells","Germinal center B cells","Plasma IgG","Plasma IgG prolif","MALT B cells","inflammation-associated B cells","Tumor-associated B cells"))
table(Idents(Bcell1))

# Set the order of cell subsets
Idents(Bcell1) <- factor(Idents(Bcell1), 
                         levels = c("Naive B cells","Follicular B cells","Germinal center B cells","Plasma IgG","Plasma IgG prolif","MALT B cells","inflammation-associated B cells","Tumor-associated B cells"))

# UMAP plot after cell annotation
pdf("03_2000_0.6_Bcell_cell_maker.pdf",width=8, height=4)
color_sub <- c("#FAFD7CFF","#82491EFF","#24325FFF","#B7E4F9FF","#FB6467FF","#526E2DFF","#E762D7FF","#E89242FF","#d3d3d3","#A6EEE6FF","#917C5DFF","#69C8ECFF")
DimPlot(Bcell1, reduction = "umap", pt.size = .1, cols = color_sub)
dev.off() 

# Save B cell Seurat object into 04Bcell directory
saveRDS(Bcell1, file.path(dir_04Bcell, "B_umap.rds"))

#==================== 2. Add response metadata (3-category Group_Curativeeffect) ====================
sample_id <- sapply(colnames(Bcell1), function(x) unlist(strsplit(x, "\\_"))[1])

# 3 categories: Resistant / Sensitive / Intermediate / Unknown
group_3 <- dplyr::case_when(
  sample_id %in% c('C100','C102','C104','C12','C50','C52','C54','C60','C62','C64','N64') ~ "Resistant",
  sample_id %in% c('C110','C122','C130','C132','C134','C20','C70','C90','C92','C94','N114','N134') ~ "Sensitive",
  sample_id %in% c('C30','C32','C40','C42','C44','C80','N44') ~ "Intermediate",
  sample_id == 'C120' ~ "Unknown",
  TRUE ~ "Sensitive"
)
Bcell1 <- AddMetaData(Bcell1, metadata = group_3, col.name = "Group_Curativeeffect")

# 2 categories: Curativeeffect (Resistant / Sensitive / Unknown)
group_2 <- dplyr::case_when(
  sample_id %in% c('C50','C60','C100','C62','C102','C54','C64','C104','C52','C80','N64') ~ "Resistant",
  sample_id == 'C120' ~ "Unknown",
  TRUE ~ "Sensitive"
)
Bcell1 <- AddMetaData(Bcell1, metadata = group_2, col.name = "Curativeeffect")
# Binary response label TWO_Curativeeffect for subsequent RO/E analysis
Bcell1 <- AddMetaData(Bcell1, metadata = group_2, col.name = "TWO_Curativeeffect")


Bcell1@meta.data$cluster_annotation <- Idents(Bcell1)
Bcell1@meta.data$subcluster <- Idents(Bcell1)

#==================== 3. Stacked bar plot: B cell proportion and absolute count per sample ====================
setwd(dir_ratio)
col_group <- c("#FAFD7CFF","#82491EFF","#24325FFF","#B7E4F9FF","#FB6467FF","#526E2DFF","#E762D7FF","#E89242FF","#d3d3d3","#A6EEE6FF","#917C5DFF","#69C8ECFF")

## 3.1 Cell proportion plot
cell.prop <- as.data.frame(prop.table(table(Idents(Bcell1), Bcell1$orig.ident)))
colnames(cell.prop) <- c("cluster","sample","proportion")
write.csv(cell.prop, "1sample_cell.prop.csv", row.names = F)

p_prop <- ggplot(cell.prop, aes(sample, proportion, fill=cluster)) +
  geom_col(stat="identity", position = "fill", width=0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_group) +
  scale_x_discrete(limits=c("C120","C20","C70","C90","C110","C130","C30","C40","C80","C50","C60","C100","C92","C122","C132","C32","C42","C52","C62","C102","C12","C94","C134","C44","C54","C64","C104","N44","N64","N94","N114","N134")) +
  theme(plot.margin=unit(rep(1,4),'cm')) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
ggsave("03B_per_sample_cell_proportion.pdf", p_prop, units="in", dpi=300, width=10, height=3, device="pdf")

## 3.2 Absolute cell count plot
cell.num2 <- as.data.frame(table(Idents(Bcell1), Bcell1$orig.ident))
colnames(cell.num2) <- c("cluster","patient","cellnum")
write.csv(cell.num2, "Bcell1_Bcell1.num.csv", row.names = F)

p_count <- ggplot(cell.num2, aes(patient, cellnum, fill=cluster)) +
  geom_col(stat="identity", width=0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_fill_manual(values = col_group) +
  scale_x_discrete(limits=c("C120","C20","C70","C90","C110","C130","C30","C40","C80","C50","C60","C100","C92","C122","C132","C32","C42","C52","C62","C102","C12","C94","C134","C44","C54","C64","C104","N44","N64","N94","N114","N134")) +
  theme(plot.margin=unit(rep(1,4),'cm')) +
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))
ggsave("03B_per_sample_cell_num.pdf", p_count, units="in", dpi=300, width=10, height=3, device="pdf")

#==================== 4. DotPlot of B cell marker genes ====================
setwd(dir_04Bcell)
markers.to.plot <- c("MS4A1","IGHD","FCER2","TCL1A","CD74","HLA-DRA","LRMP","SUGCT","SDC1","IGHG1","CD38","MKI67","TOP2A","IGHA1","JCHAIN","FTH1","TYROBP","MT1G","AGR2")
pdf("maker_Bcell1.pdf",width=10, height=4)
DotPlot(Bcell1, features = markers.to.plot, assay='RNA', group.by = 'cluster_annotation') +
  RotatedAxis()+
  theme(axis.text=element_text(size=10),
        axis.text.x = element_text(vjust = 1,size=12),
        axis.text.y = element_text(size=14))
dev.off()

#==================== 5. RORE: Chi-square residual heatmap for 3-category response ====================
setwd(dir_rore)

Bcell1_CRC_0 <- subset(Bcell1, TimePoints=="CRC_0")
Bcell1_CRC_0 <- subset(Bcell1_CRC_0, Group_Curativeeffect=="Unknown", invert=TRUE)
specificity_CRC_0 <- table(Bcell1_CRC_0$subcluster, Bcell1_CRC_0$Group_Curativeeffect)
colnames(specificity_CRC_0) <- c("Intermediate_0","Resistant_0","Sensitive_0")

Bcell1_CRC_2 <- subset(Bcell1, TimePoints=="CRC_2")
Bcell1_CRC_2 <- subset(Bcell1_CRC_2, Group_Curativeeffect=="Unknown", invert=TRUE)
specificity_CRC_2 <- table(Bcell1_CRC_2$subcluster, Bcell1_CRC_2$Group_Curativeeffect)
colnames(specificity_CRC_2) <- c("Intermediate_2","Resistant_2","Sensitive_2")

Bcell1_CRC_4 <- subset(Bcell1, TimePoints=="CRC_4")
Bcell1_CRC_4 <- subset(Bcell1_CRC_4, Group_Curativeeffect=="Unknown", invert=TRUE)
specificity_CRC_4 <- table(Bcell1_CRC_4$subcluster, Bcell1_CRC_4$Group_Curativeeffect)
colnames(specificity_CRC_4) <- c("Intermediate_4","Resistant_4","Sensitive_4")

Bcell1_CRC_Normal <- subset(Bcell1, TimePoints=="CRC_Normal")
Bcell1_CRC_Normal <- subset(Bcell1_CRC_Normal, Group_Curativeeffect=="Unknown", invert=TRUE)
specificity_CRC_Normal <- table(Bcell1_CRC_Normal$subcluster,Bcell1_CRC_Normal$TimePoints)

specificity <- cbind(specificity_CRC_0,specificity_CRC_2,specificity_CRC_4,specificity_CRC_Normal)
CMHtest(specificity)

expected <- chisq.test(specificity)$expected
result <- specificity / expected
df <- t(result)
c <- scale(df,center=T,scale=T)
all <- t(c)
write.csv(all, "Bcell1_all_R_obs_exp.csv", row.names = TRUE)
all <- read.csv("Bcell1_all_R_obs_exp.csv", row.names = 1)

pdf("01_newBcell1_Group_Curativeeffect_R_obs_exp.pdf",width=4, height=4)
annotation_col = data.frame(
  TimePoints = factor(rep(c("CRC_0","CRC_2","CRC_4","Normal"), c(3, 3, 3,1))),
  Group=factor(c(rep(c("Intermediate","Resistant","Sensitive"),3),rep("Normal",1)))
)
rownames(annotation_col) = colnames(all)
ann_colors = list(
  TimePoints = c(CRC_0 = "#008000", CRC_2 = "#0000FF", CRC_4 = "#DC143C","Normal"="#BC8F8F"),
  Group = c(Intermediate="#2F4F4F",Resistant="#B8860B",Sensitive="#008B8B",Normal="#BC8F8F")
)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pheatmap(all,
         color = c(colorRampPalette(c("#0000cc","#f5f5f5"))(length(bk)/2),colorRampPalette(c("#f5f5f5","#cc0000"))(length(bk)/2)),
         breaks=bk,legend_breaks=seq(-2,2,2),
         cluster_row = FALSE,cluster_col = FALSE,
         annotation_col = annotation_col, annotation_colors = ann_colors)
pheatmap(all,color = colorRampPalette(c("#0000cc","#f5f5f5","#cc6600","#cc0000"))(50),cluster_row = FALSE,cluster_col = FALSE,annotation_col = annotation_col, annotation_colors = ann_colors)
pheatmap(all,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cluster_row = FALSE,cluster_col = FALSE,annotation_col = annotation_col, annotation_colors = ann_colors)
dev.off()

#==================== 6. RO/E: Heatmap for binary response group (Resistant vs Sensitive) ====================
Bcell1_CRC_0 <- subset(Bcell1,TimePoints=="CRC_0")
Bcell1_CRC_0 <- subset(Bcell1_CRC_0,TWO_Curativeeffect=="Unknown",invert=TRUE)
specificity_CRC_0 <- table(Bcell1_CRC_0$subcluster, Bcell1_CRC_0$TWO_Curativeeffect)
colnames(specificity_CRC_0) <- c("Resistant_0","Sensitive_0")

Bcell1_CRC_2 <- subset(Bcell1,TimePoints=="CRC_2")
Bcell1_CRC_2 <- subset(Bcell1_CRC_2,TWO_Curativeeffect=="Unknown",invert=TRUE)
specificity_CRC_2 <- table(Bcell1_CRC_2$subcluster, Bcell1_CRC_2$TWO_Curativeeffect)
colnames(specificity_CRC_2) <- c("Resistant_2","Sensitive_2")

Bcell1_CRC_4 <- subset(Bcell1,TimePoints=="CRC_4")
Bcell1_CRC_4 <- subset(Bcell1_CRC_4,TWO_Curativeeffect=="Unknown",invert=TRUE)
specificity_CRC_4 <- table(Bcell1_CRC_4$subcluster, Bcell1_CRC_4$TWO_Curativeeffect)
colnames(specificity_CRC_4) <- c("Resistant_4","Sensitive_4")

Bcell1_CRC_Normal <- subset(Bcell1,TimePoints=="CRC_Normal")
Bcell1_CRC_Normal <- subset(Bcell1_CRC_Normal,TWO_Curativeeffect=="Unknown",invert=TRUE)
specificity_CRC_Normal <- table(Bcell1_CRC_Normal$subcluster,Bcell1_CRC_Normal$TimePoints)

specificity <- cbind(specificity_CRC_0,specificity_CRC_2,specificity_CRC_4,specificity_CRC_Normal)
CMHtest(specificity)
expected <- chisq.test(specificity)$expected
result <- specificity / expected
df <- t(result)
c <- scale(df,center=T,scale=T)
all <- t(c)
write.csv(all,"Bcell1_TWO_R_obs_exp.csv", row.names = TRUE)
all <- read.csv("Bcell1_TWO_R_obs_exp.csv",row.names = 1)

pdf("01Bcell1_TWO_Curativeeffect_R_obs_exp.pdf",width=5, height=4)
annotation_col = data.frame(
  TimePoints = factor(rep(c("CRC_0","CRC_2","CRC_4","Normal"), c(2, 2, 2,1))),
  Group=factor(c(rep(c("Resistant","Sensitive"),3),rep("Normal",1)))
)
rownames(annotation_col) = colnames(all)
ann_colors = list(
  TimePoints = c(CRC_0 = "#008000", CRC_2 = "#0000FF", CRC_4 = "#DC143C","Normal"="#BC8F8F"),
  Group = c("Resistant"="#B8860B","Sensitive"="#2F4F4F","Normal"="#BC8F8F")
)
pheatmap(all,color = colorRampPalette(c("#0000cc","#f5f5f5","#cc6600","#cc0000"))(50),cluster_row = FALSE,cluster_col = FALSE,annotation_col = annotation_col, annotation_colors = ann_colors)
pheatmap(all,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),cluster_row = FALSE,cluster_col = FALSE,annotation_col = annotation_col, annotation_colors = ann_colors)
dev.off()