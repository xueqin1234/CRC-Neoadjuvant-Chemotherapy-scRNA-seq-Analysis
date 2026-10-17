#==================== Load packages ====================
library(tidyverse)
library(DT)

#==================== Global Functions ====================
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                options = list(dom = 'Blfrtip',
                               buttons = c('csv', 'excel', 'pdf'),
                               lengthMenu = list(c(10,25,50,-1), c(10,25,50,"All"))))
}

# Optimized order_sequence: vectorized instead of slow for-loop
order_sequence <- function(df){
  # Only need swap gene_a / gene_b when receptor_a=TRUE & receptor_b=FALSE
  swap_idx <- df$receptor_a == "True" & df$receptor_b == "False"
  if(sum(swap_idx) > 0){
    sub_swap <- df[swap_idx, ]
    old_colnames <- colnames(sub_swap)
    # swap gene_a and gene_b
    sub_swap <- sub_swap %>% mutate(gene_a_tmp = gene_a, gene_a = gene_b, gene_b = gene_a_tmp) %>% select(-gene_a_tmp)
    # swap the cell pair columns (col12 onwards: cellA|cellB --> cellB|cellA)
    cell_cols <- old_colnames[12:length(old_colnames)]
    new_cell_names <- map_chr(cell_cols, ~{
      sp <- strsplit(.x, "\\|")[[1]]
      paste(sp[2], sp[1], sep = "|")
    })
    colnames(sub_swap)[12:length(old_colnames)] <- new_cell_names
    df[swap_idx,] <- sub_swap
  }
  return(df)
}

#==================== CPDB pipeline function for single timepoint ====================
cpdb_analysis <- function(out_dir, time_tag){
  # read cpdb output
  mymeans <- read.table(file.path(out_dir, "means.txt"), header = T, sep = "\t", check.names=F)
  mypvals <- read.table(file.path(out_dir, "pvalues.txt"), header = T, sep = "\t", check.names=F)
  
  # Filter L-R pairs: one receptor, one ligand
  df <- subset(mymeans, (receptor_a == 'True' & receptor_b == 'False') | (receptor_a == 'False' & receptor_b == 'True'))
  # Remove rows with NA/empty
  df <- df %>% mutate(na_count = rowSums(is.na(.) | . == "")) %>% filter(na_count == 0) %>% select(-na_count)
  df <- df[!duplicated(df$interacting_pair), ]
  mypvals <- mypvals[!duplicated(mypvals$interacting_pair), ]
  
  # reorder geneA/geneB
  means_order <- order_sequence(df) %>% unite(Pairs, gene_a, gene_b)
  pvals_order <- order_sequence(mypvals) %>% unite(Pairs, gene_a, gene_b)
  
  # subset columns
  means_sub <- means_order[, c('Pairs', colnames(mymeans)[-c(1:11)])]
  pvals_sub <- pvals_order[, c('Pairs', colnames(mymeans)[-c(1:11)])]
  
  # longer format (replace gather)
  means_gather <- means_sub %>% pivot_longer(cols=-Pairs, names_to="celltype", values_to="mean_expression")
  pvals_gather <- pvals_sub %>% pivot_longer(cols=-Pairs, names_to="celltype", values_to="pval")
  
  mean_pval <- left_join(means_gather, pvals_gather, by = c("Pairs", "celltype"))
  
  # Filter significant pairs: at least one cell pair pval <=0.05
  a <- mean_pval %>% select(Pairs, celltype, pval) %>% pivot_wider(names_from = celltype, values_from = pval)
  sig_pairs <- a[rowSums(a[, -1] <= 0.05, na.rm = TRUE) !=0, ]
  
  mean_pval_sub <- mean_pval %>% filter(Pairs %in% sig_pairs$Pairs)
  mean_pval_sub$Timepoint <- time_tag
  return(mean_pval_sub)
}

#==================== Path Configuration ====================
root_dir <- "/hanlab/xueqinliu/test"
cpdb_base <- file.path(root_dir, "cellphoneDB_to_CRC/cellphoneDB_laoban/Intermediate")

# output folders for each timepoint
cpdb_out0 <- file.path(cpdb_base, "CRC0/out")
cpdb_out2 <- file.path(cpdb_base, "CRC2/out")
cpdb_out4 <- file.path(cpdb_base, "CRC4/out")

# define your interested cell pairs ONCE
cell.pair <- c("malignant|neutrophil","malignant|FAP+ fibro","malignant|Tip-like ECs","malignant|Pericyte",
               "neutrophil|malignant","FAP+ fibro|malignant","Tip-like ECs|malignant","Pericyte|malignant",
               "neutrophil|FAP+ fibro","neutrophil|Tip-like ECs","neutrophil|Pericyte",
               "FAP+ fibro|neutrophil","Tip-like ECs|neutrophil","Pericyte|neutrophil")

#==================== Run CPDB analysis for 3 time points ====================
mean_pval_sub_CRC0 <- cpdb_analysis(cpdb_out0, "CRC_0")
mean_pval_sub_CRC2 <- cpdb_analysis(cpdb_out2, "CRC_2")
mean_pval_sub_CRC4 <- cpdb_analysis(cpdb_out4, "CRC_4")

# Merge all timepoints
finally <- bind_rows(mean_pval_sub_CRC0, mean_pval_sub_CRC2, mean_pval_sub_CRC4)

# Filter expression threshold
data <- finally %>% filter(mean_expression > 0.2)
data1 <- data %>% filter(celltype %in% cell.pair)

# Pivot pvalue table and export csv
filter_table <- data1 %>% select(Pairs, celltype, mean_expression, Timepoint, pval) %>%
  pivot_wider(names_from = Timepoint, values_from = pval)

# Fill NA pvalue to 1
filter_table <- filter_table %>% mutate(across(where(is.numeric), ~replace_na(.x, 1)))
write.csv(filter_table, file.path(root_dir, "myeloid_stromal.csv"), row.names = F)

data2 <- data1 %>% filter(Pairs %in% filter_table$Pairs & celltype %in% filter_table$celltype)

# Set fixed order for x-axis cell pairs
cell_pair_levels <- cell.pair
data2$celltype <- factor(data2$celltype, levels = cell_pair_levels)

#==================== Plot & save PDF ====================
pdf(file.path(root_dir, "02candidate_Intermediate.pdf"), width=10, height=9)

ggplot(data2,aes(celltype,Pairs))+
  geom_point(aes(color=log2(mean_expression+1e-3), size=-log10(pval+1e-3))) +
  facet_wrap(~Timepoint, ncol = 3)+
  scale_size_continuous(range = c(0,2))+
  scale_color_gradient2(high="red",mid = "yellow",low ="darkblue")+
  theme_bw()+
  theme(
    axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, face="italic",size=9),
    axis.text.y = element_text(size=8)
  )

dev.off()

# Optional: interactive table
# create_dt(data2)