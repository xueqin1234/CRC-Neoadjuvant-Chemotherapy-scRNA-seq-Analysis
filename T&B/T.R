#==================== 加载包 ====================
library(Seurat)
library(harmony)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(vcdExtra)
library(vcd)
library(dplyr) # case_when依赖，原代码缺失

#==================== 【核心路径配置，只需要在这里修改】 ====================
root_dir <- "/hanlab/xueqinliu/test"
# 子文件夹
dir_03tcell <- file.path(root_dir, "03tcell")
dir_ratio   <- file.path(dir_03tcell, "03ratio")
dir_rore    <- file.path(dir_03tcell, "02RORE")

# 自动创建文件夹，不存在则新建
if (!dir.exists(root_dir)) dir.create(root_dir, recursive = TRUE)
if (!dir.exists(dir_03tcell)) dir.create(dir_03tcell, recursive = TRUE)
if (!dir.exists(dir_ratio)) dir.create(dir_ratio, recursive = TRUE)
if (!dir.exists(dir_rore)) dir.create(dir_rore, recursive = TRUE)

# 原始数据存放目录（你原始rds还在这个路径）
data_raw_dir <- "/hanlab/xueqin/03zhoulian/data"

#==================== 1. 读取整合后CRC单细胞，提取T细胞子集并重新分析 ====================
setwd(data_raw_dir)
allcell <- readRDS("03CRC_ALL_single_harmony.rds")

# 重命名cluster -> 大细胞类型
allcell <- RenameIdents(allcell, 
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

# 剔除Basal、Acinar；剔除C114样本
allcell_filter <- subset(allcell, idents = c("Basal", "Acinar cells"), invert = TRUE)
allcell_filter <- subset(allcell_filter, orig.ident == "C114", invert = TRUE)

# 设置细胞类型顺序
Idents(allcell_filter) <- factor(Idents(allcell_filter), 
                                 levels = c("Epithelial cells","Stromal cells","T cells","B cells","Myeloid","Mast"))
table(Idents(allcell_filter))

# 提取T细胞子集
tcell_raw <- subset(allcell_filter, idents = "T cells")

# 重建seurat对象 + 去除线粒体基因
tcell <- CreateSeuratObject(counts = tcell_raw@assays$RNA@counts, meta.data = tcell_raw@meta.data)
mt.index <- grep("^MT-", rownames(tcell), value = F)
tcell <- tcell[-mt.index, ]

# 标准单细胞流程
tcell <- NormalizeData(tcell, verbose = FALSE)
tcell <- FindVariableFeatures(tcell, selection.method = 'vst', nfeatures = 2000)
hvgs <- VariableFeatures(tcell)
tcell <- ScaleData(tcell, features = hvgs, verbose = FALSE)
tcell <- RunPCA(tcell, features = hvgs, npcs = 80, verbose = FALSE)

# 可选：判断PCA维度
# ElbowPlot(tcell, ndims = 50)
# DimHeatmap(tcell, dims = 1:50, balanced = TRUE)
tcell <- RunUMAP(tcell, reduction = "pca", dims = 1:40, verbose = FALSE)
tcell <- FindNeighbors(tcell, reduction = "pca", dims = 1:40, verbose = FALSE)
tcell <- FindClusters(tcell, resolution = 0.5, verbose = FALSE)

# 绘图：原始聚类UMAP（输出到03tcell）
setwd(dir_03tcell)
pdf("01_Tcell_raw_cluster_umap.pdf", width=7, height=5)
DimPlot(tcell, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
dev.off()

#==================== 2. T细胞亚群重命名，剔除混杂非T细胞（B/mono） ====================
tcell_anno <- RenameIdents(tcell, 
                           `0` = "Naïve", 
                           `1` = "CD8+ Cytotoxic", 
                           `2` = "Treg cells", 
                           `3` = "Effector memory cells", 
                           `4` = "MALT B cells", 
                           `5` = "Th 17 cells", 
                           `6` = "CD8+ T cells", 
                           `7` = "HSP+ T cells", 
                           `8` = "CD8+ Cytotoxic", 
                           `9` = "CD8+ T cells", 
                           `10` = "CD8+ MAIT cells", 
                           `11` = "CD4+ T cells", 
                           `12` = "Naïve", 
                           `13` = "NK cells", 
                           `14` = "MALT B cells",
                           `15` = "CD8+ Exhausted", 
                           `16` = "T follicular helper",
                           `17` = "gamma delta T", 
                           `18` = "MALT B cells",
                           `19` = "IFI T cells", 
                           `20` = "ILC3 T",
                           `21` = "mono")

# 去掉MALT B / mono污染细胞
tcell_filtered <- subset(tcell_anno, 
                         idents = c("Naïve","Effector memory cells","CD4+ T cells","Treg cells","Th 17 cells",
                                    "T follicular helper","CD8+ T cells","CD8+ Cytotoxic","CD8+ Exhausted",
                                    "CD8+ MAIT cells","NK cells","ILC3 T","gamma delta T","IFI T cells","HSP+ T cells"))

# 设置亚群顺序
subcluster_level <- c("Naïve","Effector memory cells","CD4+ T cells","Treg cells","Th 17 cells",
                      "T follicular helper","CD8+ T cells","CD8+ Cytotoxic","CD8+ Exhausted",
                      "CD8+ MAIT cells","NK cells","ILC3 T","gamma delta T","IFI T cells","HSP+ T cells")
Idents(tcell_filtered) <- factor(Idents(tcell_filtered), levels = subcluster_level)
table(Idents(tcell_filtered))

# 亚群注释UMAP图
color_sub <- c("#336633","#6b8e23","#adff2f","#3cb371","#0000ff","#0066ff","#008b8b","#9999ff",
               "#00ced1","#6495ed","#f4a460","#f0e68c","#d2691e","#ff6347","#b22222")

pdf("02sec_sub_Tcell_cellnames.pdf", width=7, height=5)
DimPlot(tcell_filtered, reduction = "umap", pt.size = 0.5, cols = color_sub)
dev.off()

# 保存T细胞对象，放到test根目录
saveRDS(tcell_filtered, file.path(root_dir, "tcell_umap.rds"))

# 提取样本名前缀
sample_id <- sapply(colnames(tcell_filtered), function(x) unlist(strsplit(x, "\\_"))[1])

## 分组1：Group_Curativeeffect（三分类：Sensitive / Resistant / Intermediate / Unknown）
group1 <- dplyr::case_when(
  sample_id %in% c('C100','C102','C104','C12','C50','C52','C54','C60','C62','C64','N64') ~ "Resistant",
  sample_id %in% c('C110','C122','C130','C132','C134','C20','C70','C90','C92','C94','N114','N134') ~ "Sensitive",
  sample_id %in% c('C30','C32','C40','C42','C44','C80','N44') ~ "Intermediate",
  sample_id == 'C120' ~ "Unknown",
  TRUE ~ "Sensitive"
)
tcell_filtered <- AddMetaData(tcell_filtered, metadata = group1, col.name = 'Group_Curativeeffect')

## 分组2：Curativeeffect（二分类，你原来的第二套规则）
group2 <- dplyr::case_when(
  sample_id %in% c('C50','C60','C100','C62','C102','C54','C64','C104','C52','C80','N64') ~ "Resistant",
  sample_id == 'C120' ~ "Unknown",
  TRUE ~ "Sensitive"
)
tcell_filtered <- AddMetaData(tcell_filtered, metadata = group2, col.name = 'Curativeeffect')

tcell_filtered@meta.data$subcluster <- Idents(tcell_filtered)
tcell_filtered@meta.data$cluster_annotation <- Idents(tcell_filtered)

#==================== 4. 样本水平细胞比例&细胞数 堆叠柱状图 ====================
setwd(dir_ratio)

# 4.1 每个样本内T亚群占比
cell_prop_df <- as.data.frame(prop.table(table(Idents(tcell_filtered), tcell_filtered$orig.ident)))
colnames(cell_prop_df) <- c("cluster","sample","proportion")
write.csv(cell_prop_df, "1sample_cell.prop.csv", row.names = F)

# 堆叠百分比柱状图
p_prop <- ggplot(cell_prop_df, aes(x = sample, y = proportion, fill = cluster)) +
  geom_col(position = "fill", width = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(rep(1,4), 'cm'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = color_sub) +
  scale_x_discrete(limits=c("C120","C20","C70","C90","C110","C130","C30","C40","C80","C50","C60","C100",
                            "C92","C122","C132","C32","C42","C52","C62","C102","C12","C94","C134","C44",
                            "C54","C64","C104","N44","N64","N94","N114","N134"))

ggsave("03T_per_sample_cell_proportion.pdf", p_prop, units="in", dpi=300, width=10, height=3, device="pdf")

# 4.2 每个样本各亚群细胞绝对数量
cell_count_df <- as.data.frame(table(Idents(tcell_filtered), tcell_filtered$orig.ident))
colnames(cell_count_df) <- c("cluster","patient","cellnum")
write.csv(cell_count_df, "tcell_filtered_tcell_filtered.num.csv", row.names = F)

p_count <- ggplot(cell_count_df, aes(x = patient, y = cellnum, fill = cluster)) +
  geom_col(width = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(rep(1,4), 'cm'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = color_sub) +
  scale_x_discrete(limits=c("C120","C20","C70","C90","C110","C130","C30","C40","C80","C50","C60","C100",
                            "C92","C122","C132","C32","C42","C52","C62","C102","C12","C94","C134","C44",
                            "C54","C64","C104","N44","N64","N94","N114","N134"))

ggsave("03T_per_sample_cell_num.pdf", p_count, units="in", dpi=300, width=10, height=3, device="pdf")


#==================== 6. DotPlot 标记基因点图 ====================
markers.to.plot <- c("CCR7","TCF7","SELL","CD44","PRDM1","IL7R","CD4","FOXP3","IL2RA","IKZF2","IL17A","IL17F","RORA",
                     "CXCL13","MAF","CXCR5","CD8A","CD8B","IFNG","GZMA","NKG7","HAVCR2","LAG3","CTLA4","ZBTB16",
                     "SLC4A10","FCGR3A","KLRF1","KLRD1","KIT","IL1R1","IL23R","TRDC","TRGC2","TRGC1","IFI6","IFIT3",
                     "IFI44L","HSPA1A","HSPH1")

pdf("02maker_tcell_filtered.pdf", width=12, height=5)
DotPlot(tcell_filtered, features = markers.to.plot, assay='RNA', group.by = 'cluster_annotation') +
  RotatedAxis() +
  theme(axis.text = element_text(size=7),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(face="italic",size=11))
dev.off()

#==================== 7. 不同时间点疗效分组列联表，卡方残差热图 ====================
setwd(dir_rore)

# 分时间点提取，剔除Unknown
T_CRC_0 <- subset(tcell_filtered, TimePoints=="CRC_0")
T_CRC_0 <- subset(T_CRC_0, Group_Curativeeffect=="Unknown", invert=TRUE)
tab0 <- table(T_CRC_0$subcluster, T_CRC_0$Group_Curativeeffect)
colnames(tab0) <- c("Intermediate_0","Resistant_0","Sensitive_0")

T_CRC_2 <- subset(tcell_filtered, TimePoints=="CRC_2")
T_CRC_2 <- subset(T_CRC_2, Group_Curativeeffect=="Unknown", invert=TRUE)
tab2 <- table(T_CRC_2$subcluster, T_CRC_2$Group_Curativeeffect)
colnames(tab2) <- c("Intermediate_2","Resistant_2","Sensitive_2")

T_CRC_4 <- subset(tcell_filtered, TimePoints=="CRC_4")
T_CRC_4 <- subset(T_CRC_4, Group_Curativeeffect=="Unknown", invert=TRUE)
tab4 <- table(T_CRC_4$subcluster, T_CRC_4$Group_Curativeeffect)
colnames(tab4) <- c("Intermediate_4","Resistant_4","Sensitive_4")

T_CRC_Normal <- subset(tcell_filtered, TimePoints=="CRC_Normal")
T_CRC_Normal <- subset(T_CRC_Normal, Group_Curativeeffect=="Unknown", invert=TRUE)
tabN <- table(T_CRC_Normal$subcluster, T_CRC_Normal$TimePoints)

# 合并矩阵
spec_mat <- cbind(tab0, tab2, tab4, tabN)

# CMH检验
CMHtest(spec_mat)

# 卡方：观测值 / 期望值
chi_res <- chisq.test(spec_mat)
exp_mat <- chi_res$expected
obs_over_exp <- spec_mat / exp_mat

# Z-score标准化
df_t <- t(obs_over_exp)
df_scale <- scale(df_t, center = T, scale = T)
df_final <- t(df_scale)

write.csv(df_final, "Tcell_all_R_obs_exp.csv")

# 热图注释
annotation_col = data.frame(
  TimePoints = factor(rep(c("CRC_0","CRC_2","CRC_4","Normal"), c(3, 3, 3,1))),
  Group = factor(c(rep(c("Intermediate","Resistant","Sensitive"),3), rep("Normal",1)))
)
rownames(annotation_col) = colnames(df_final)

ann_colors = list( 
  TimePoints = c(CRC_0 = "#008000", CRC_2 = "#0000FF", CRC_4 = "#DC143C","Normal"="#BC8F8F"),
  Group = c(Intermediate="#2F4F4F",Resistant="#B8860B",Sensitive="#008B8B",Normal="#BC8F8F")
)

# 绘图
pdf("01_newTcell_Group_Curativeeffect_R_obs_exp.pdf", width=5, height=4)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
pheatmap(df_final,
         color = c(colorRampPalette(c("#0000cc","#f5f5f5"))(length(bk)/2),
                   colorRampPalette(c("#f5f5f5","#cc0000"))(length(bk)/2)),
         breaks=bk,
         legend_breaks=seq(-2,2,2),
         cluster_row = FALSE,
         cluster_col = FALSE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors)

# 备用配色（三色蓝白红，推荐用这个，更符合单细胞习惯）
pheatmap(df_final,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         cluster_row = FALSE,
         cluster_col = FALSE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors)
dev.off()