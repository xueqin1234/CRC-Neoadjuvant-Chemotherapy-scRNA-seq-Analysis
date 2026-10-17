#=========================================================================
# NicheNet analysis script: Resistant CRC0, neutrophil( sender ) -> FAP+ fibro( receiver )
# Working root: /hanlab/xueqinliu/test
# Optimized: unified file paths, case_when instead of nested ifelse, reusable function, clean redundant comments and duplicate code
#=========================================================================
library(Seurat)
library(dplyr)
library(tidyverse)
library(nichenetr)
library(circlize)
library(cowplot)
library(ggpubr)
library(RColorBrewer)

#==================== Path configuration section (modify paths only here) ====================
root_dir       <- "./test"
raw_data_dir   <- "./test/data"
external_rds   <- "./test/00data"
ref_dir        <- "./ref"
nichenet_out   <- file.path(root_dir, "02nichenet/res")
dir.create(nichenet_out, recursive = TRUE, showWarnings = FALSE)
#==================== Load main Seurat object and add response metadata ====================
allcel <- readRDS(file.path(external_rds, "all_cell_ano.rds"))
# Extract Resistant group and split by time points
Resistant <- subset(allcel, Group_Curativeeffect == "Resistant")
Resistant_crc0 <- subset(Resistant, TimePoints == "CRC_0")
Resistant_crc2 <- subset(Resistant, TimePoints == "CRC_2")
Resistant_crc4 <- subset(Resistant, TimePoints == "CRC_4")

#==================== Core NicheNet reusable function (sender -> receiver) ====================
run_nichenet <- function(seu_obj, sender_celltype, receiver_celltype, ligand_target_matrix, lr_network, out_dir){
  # 1. Get expressed genes of sender and receiver cells
  expressed_genes_receiver <- get_expressed_genes(receiver_celltype, seu_obj, pct = 0.20)
  list_expressed_genes_sender <- lapply(unique(sender_celltype), function(x) get_expressed_genes(x, seu_obj, pct = 0.10))
  expressed_genes_sender <- unlist(list_expressed_genes_sender) %>% unique()

  # 2. DEG analysis for Receiver(FAP+ fibro): Resistant vs Sensitive
  seurat_obj_receiver <- subset(seu_obj, idents == receiver_celltype)
  Idents(seurat_obj_receiver) <- seurat_obj_receiver$Group_Curativeeffect
  DE_table_receiver <- FindMarkers(
    object = seurat_obj_receiver,
    ident.1 = "Resistant",
    ident.2 = "Sensitive",
    min.pct = 0.10,
    logfc.threshold = 0.25,
    only.pos = TRUE
  )
  DE_table_receiver$gene <- rownames(DE_table_receiver)
  geneset_oi <- DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.35) %>% pull(gene)
  geneset_oi <- intersect(geneset_oi, rownames(ligand_target_matrix))

  # Background expressed genes
  background_expressed_genes <- intersect(expressed_genes_receiver, rownames(ligand_target_matrix))

  # Filter potential ligand-receptor pairs
  ligands <- lr_network %>% pull(from) %>% unique()
  expressed_ligands <- intersect(ligands, expressed_genes_sender)
  receptors <- lr_network %>% pull(to) %>% unique()
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)
  potential_ligands <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

  # Predict ligand activities
  ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi,
    background_expressed_genes = background_expressed_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = potential_ligands
  )
  best_upstream_ligands <- ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)

  # Ligand-target gene network
  active_ligand_target_links_df <- lapply(
    best_upstream_ligands,
    get_weighted_ligand_target_links,
    geneset = geneset_oi,
    ligand_target_matrix = ligand_target_matrix,
    n = 250
  ) %>% bind_rows() %>% na.omit()

  active_ligand_target_links <- prepare_ligand_target_visualization(
    ligand_target_df = active_ligand_target_links_df,
    ligand_target_matrix = ligand_target_matrix,
    cutoff = 0.25
  )

  order_ligands <- rev(intersect(best_upstream_ligands, colnames(active_ligand_target_links)))
  order_targets <- unique(active_ligand_target_links_df$target)
  vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])

  write.csv(vis_ligand_target, file.path(out_dir, "neutrophil_ligands_FAP_target1.csv"))

  # Target gene heatmap
  p_ligand_target_network <- vis_ligand_target %>%
    make_heatmap_ggplot("Prioritized neutrophil ligands","genes in FAP+ fibro cells",
                        color = "purple",legend_position = "top", x_axis_position = "top",
                        legend_title = "Regulatory potential") +
    scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) +
    theme(axis.text.x = element_text(face = "italic"))

  # Pearson correlation heatmap
  ligand_pearson_matrix <- ligand_activities %>% select(pearson) %>%
    as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
  vis_ligand_pearson <- ligand_pearson_matrix[order_ligands, , drop=F] %>% magrittr::set_colnames("Pearson")
  p_ligand_pearson <- vis_ligand_pearson %>%
    make_heatmap_ggplot("Prioritized neutrophil ligands","genes in FAP+ fibro cells",
                        color = "darkorange",legend_position = "top", x_axis_position = "top",
                        legend_title = "Pearson correlation coefficient\ntarget gene prediction ability")

  # Dotplot: ligand expression in sender cells
  seurat_sender <- subset(seu_obj, tumor %in% sender_celltype)
  seurat_sender <- subset(seurat_sender, !Group_Curativeeffect %in% c("Unknown","Intermediate"))
  seurat_sender <- subset(seurat_sender, TimePoints=="CRC_0")
  Idents(seurat_sender) <- seurat_sender$Group_Curativeeffect

  rotated_dotplot <- DotPlot(seurat_sender, features = rev(best_upstream_ligands), cols = "RdYlBu") +
    RotatedAxis() + coord_flip()

  # Ligand-receptor interaction heatmap
  lr_network_top <- lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
  best_upstream_receptors <- lr_network_top %>% pull(to) %>% unique()

  weighted_networks <- readRDS(file.path(ref_dir, "weighted_networks.rds"))
  lr_network_top_df <- weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
  lr_network_top_df <- lr_network_top_df %>% pivot_wider(names_from = from, values_from = weight, fill = 0)
  lr_network_top_matrix <- lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

  dist_receptors <- dist(lr_network_top_matrix, method = "binary")
  hclust_receptors <- hclust(dist_receptors, method = "ward.D2")
  order_receptors <- hclust_receptors$labels[hclust_receptors$order]

  dist_ligands <- dist(t(lr_network_top_matrix), method = "binary")
  hclust_ligands <- hclust(dist_ligands, method = "ward.D2")
  order_ligands_receptor <- hclust_ligands$labels[hclust_ligands$order]

  vis_ligand_receptor_network <- lr_network_top_matrix[order_receptors, order_ligands_receptor]
  p_ligand_receptor_network <- t(vis_ligand_receptor_network) %>%
    make_heatmap_ggplot("Prioritized neutrophil ligands","genes in FAP+ fibro cells",
                        color = "mediumvioletred", x_axis_position = "top",
                        legend_title = "Prior interaction potential")

  # Combine all panels into one figure
  figures_without_legend <- cowplot::plot_grid(
    p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
    rotated_dotplot + theme(legend.position = "none", axis.ticks = element_blank(),
                            axis.title.x = element_text(size = 12),
                            axis.text.y = element_text(face = "italic", size = 9),
                            axis.text.x = element_text(size = 9,  angle = 90,hjust = 0)) +
      ylab("Expression in Sender") + xlab("") + scale_y_discrete(position = "right"),
    p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
    align = "hv", nrow = 1,
    rel_widths = c(ncol(vis_ligand_pearson)+5,ncol(vis_ligand_pearson)+7, ncol(vis_ligand_target)))

  legends <- cowplot::plot_grid(
    ggpubr::as_ggplot(get_legend(p_ligand_pearson)),
    ggpubr::as_ggplot(get_legend(rotated_dotplot)),
    ggpubr::as_ggplot(get_legend(p_ligand_target_network)),
    nrow = 1, align = "h", rel_widths = c(1.5, 1, 1, 1))

  combined_plot <- cowplot::plot_grid(figures_without_legend, legends, rel_heights = c(10,5), nrow = 2, align = "hv")

  # Export PDF figures
  ggsave(file.path(out_dir, "01neutrophil_ligands_FAP_target-genes.pdf"), p_ligand_target_network, width = 12,height = 12)
  ggsave(file.path(out_dir, "02FAP_neutrophil_ligands_Receptors.pdf"), p_ligand_receptor_network, width = 12,height = 15)
  ggsave(file.path(out_dir, "04-3all.pdf"), combined_plot, width = 15,height = 12)

  return(list(
    ligand_activities = ligand_activities,
    best_ligands = best_upstream_ligands,
    p_target = p_ligand_target_network,
    p_pearson = p_ligand_pearson,
    p_lr = p_ligand_receptor_network,
    p_dot = rotated_dotplot,
    combined = combined_plot
  ))
}

#==================== Load NicheNet reference database and run analysis ====================
ligand_target_matrix <- readRDS(file.path(ref_dir, "ligand_target_matrix.rds"))
lr_network <- readRDS(file.path(ref_dir, "lr_network.rds"))

# Run analysis: sender=neutrophil, receiver=FAP+ fibro, Resistant_CRC0
nichenet_result <- run_nichenet(
  seu_obj = Resistant_crc0,
  sender_celltype = "neutrophil",
  receiver_celltype = "FAP+ fibro",
  ligand_target_matrix = ligand_target_matrix,
  lr_network = lr_network,
  out_dir = nichenet_out
)