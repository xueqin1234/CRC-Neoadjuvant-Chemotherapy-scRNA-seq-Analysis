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

# ===================== 路径配置，全部不用setwd =====================
root_dir <- "./test"
path_all_epi <- file.path(root_dir, "00data", "all_epi_nomalignant_malignant.rds")
out_nmf_main <- file.path(root_dir, "04epi", "01manlignant_NMF")
out_nmf_rank8 <- file.path(root_dir, "04epi", "01manlignant_NMF", "nmf_rank_8")
save_rds_path <- file.path(root_dir, "00data", "01_rank8_malignant.nmf.rds")

# 读取上皮对象，提取恶性细胞
all_epi <- readRDS(path_all_epi)
malignant <- subset(all_epi, tomur == "malignant")
head(malignant@meta.data)

# 重新构建Seurat对象，去除线粒体基因，标准化
malignant <- CreateSeuratObject(counts = malignant@assays$RNA@counts, meta.data = malignant@meta.data)
malignant[["percent.mt"]] <- PercentageFeatureSet(malignant, pattern = "^MT-")
mt.index <- grep(pattern = "^MT-", x = rownames(malignant@assays$RNA@counts), value = FALSE)
malignant <- malignant[-mt.index, ]

malignant <- NormalizeData(malignant, verbose = FALSE)
GetAssay(malignant, assay = "RNA")

# PCA流程
malignant <- FindVariableFeatures(malignant, selection.method = 'vst')
var_genes <- VariableFeatures(malignant)
malignant <- ScaleData(malignant, features = var_genes, do.center = FALSE)
malignant <- RunPCA(object = malignant, pc.genes = var_genes, npcs = 80)

# 肘部图
pdf(file.path(out_nmf_main, "ElbowPlot.pdf"), width = 6, height = 4)
print(ElbowPlot(malignant, ndims = 50))
dev.off()

# 聚类 + UMAP based on PCA
malignant <- FindNeighbors(malignant, dims = 1:35)
res.used <- 0.23
malignant <- FindClusters(object = malignant, verbose = TRUE, resolution = res.used)
table(malignant$seurat_clusters)

set.seed(123)
malignant <- RunUMAP(malignant, reduction = "pca", dims = 1:35)

# NMF输入矩阵
vm <- malignant@assays$RNA@scale.data
dim(vm)

i0 <- which(colSums(vm) == 0)
i_na <- which(colSums(is.na(vm)) > 0)
vm_filtered <- vm[, -c(i0, i_na)]
dim(vm_filtered)

# NMF rank=8
set.seed(123)
res <- nmf(vm, 8, nrun = 5, method = "snmf/r", seed = 'nndsvd')
save(res, file = file.path(out_nmf_main, "res8.RData"))


# 加载NMF结果，把NMF嵌入放到reduction
load(file.path(out_nmf_main, "res8.RData"))
malignant@reductions$nmf <- malignant@reductions$pca
malignant@reductions$nmf@cell.embeddings <- t(coef(res))
malignant@reductions$nmf@feature.loadings <- basis(res)

# NMF空间UMAP + 聚类
set.seed(219)
malignant.nmf <- RunUMAP(malignant, reduction = 'nmf', dims = 1:8) %>%
  FindNeighbors(reduction = 'nmf', dims = 1:8) %>%
  FindClusters()

# NMF最大系数作为cluster
malignant.nmf$cluster <- apply(NMF::coefficients(res)[1:8, ], 2, which.max)

# UMAP可视化pdf
pdf(file.path(out_nmf_rank8, "03malignant.nmf.pdf"))
print(DimPlot(malignant.nmf, group.by = "cluster") + scale_color_nejm())
print(DimPlot(malignant.nmf, group.by = "subcluster",
              cols = c("#f4a460","#ffd700","#90ee90","#48d1cc","#6495ed","#4169e1","#9370db","#d2691e","#db7093","#dcdcdc")))
dev.off()

# 提取NMF top300特征基因
f <- extractFeatures(res, 300L)
f <- lapply(f, function(x) rownames(res)[x])
f <- do.call("rbind", f)
DT::datatable(t(f))

write.csv(t(f), file.path(out_nmf_rank8, "top300.csv"), row.names = FALSE)
write.csv(table(malignant.nmf$orig.ident, malignant.nmf$cluster),
          file.path(out_nmf_rank8, "rank8_number.csv"))

# 分面UMAP图
pdf(file.path(out_nmf_rank8, "04malignantnmf_groupby_subcluster.pdf"))
print(DimPlot(malignant.nmf,
              group.by = "subcluster",
              cols = c("#f4a460","#ffd700","#90ee90","#48d1cc","#6495ed","#4169e1","#9370db","#d2691e","#db7093","#dcdcdc"),
              split.by = "Group_Curativeeffect", ncol = 3) +
        theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")))
dev.off()

# 保存nmf对象
saveRDS(malignant.nmf, file = save_rds_path)

# ---------------------- 批量绘图函数：FeaturePlot ----------------------
plot_feature_gene <- function(gene_name, obj, out_dir, pt_size=0.8, w=5, h=4){
  p <- FeaturePlot(
    obj,
    features = gene_name,
    pt.size = pt_size
  ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13),
      legend.position = "right"
    )
  ggsave(file.path(out_dir, paste0(gene_name,"_FeaturePlot_pt1.2.pdf")),
         p, width = w, height = h, dpi = 300)
  return(p)
}

feature_gene_list <- c("CXCL3","CXCL2","CXCL1","CA9","TDO2")
for(g in feature_gene_list){
  plot_feature_gene(g, malignant.nmf, out_nmf_rank8)
}

# ---------------------- 批量绘图函数：VlnPlot ----------------------
plot_vln_gene <- function(gene_name, obj, group_var, out_dir, pt_size=0.15, w=6, h=4){
  p <- VlnPlot(
    obj,
    group.by = group_var,
    features = gene_name,
    pt.size = pt_size
  ) +
    scale_fill_nejm() +
    scale_color_nejm() +
    labs(
      x = group_var,
      y = paste0(gene_name, " expression"),
      title = gene_name
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size =14, face = "bold"),
      axis.text.x = element_text(size =11, angle =45, hjust =1),
      axis.text.y = element_text(size =11),
      axis.title = element_text(size =12),
      legend.position = "none",
      panel.grid = element_blank()
    )
  if(group_var == "cluster"){
    p <- p + geom_boxplot(width=0.1, fill="white", outlier.size=0.2)
  }
  ggsave(file.path(out_dir, paste0(gene_name,"_",group_var,".pdf")), p, width = w, height = h, dpi =300)
  return(p)
}

plot_vln_gene("SOX9", malignant.nmf, "cluster", out_nmf_rank8)
plot_vln_gene("TDO2", malignant.nmf, "subcluster", out_nmf_rank8)