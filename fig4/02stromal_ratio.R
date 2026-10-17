# Load required packages
library(Seurat)
library(harmony)
library(ggplot2)
library(ggsci)

# Define paths
obj_path <- "./test/00data/stromal_umap_new.rds"
out_dir  <- "./test/02stromal/02ratio"

# Load stromal Seurat object
stromal <- readRDS(file = obj_path)

##===================== 1. Stacked bar plot: stromal subcluster proportion per sample =====================
# Print cell count and proportion of each stromal subcluster across all cells
table(Idents(stromal))
prop.table(table(Idents(stromal)))

# Create contingency table: subcluster vs sample
tab_count <- table(Idents(stromal), stromal@meta.data$orig.ident)
# Normalize by sample (each sample sums to 1)
cell_proportion <- prop.table(tab_count, margin = 2)
# Convert table to long format dataframe for ggplot
cell.prop <- as.data.frame(cell_proportion)
colnames(cell.prop) <- c("cluster","sample","proportion")

# Export proportion table to csv
write.csv(cell.prop, file.path(out_dir, "1sample_cell.prop.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(cell_proportion), file.path(out_dir, "1sample_cell.prop2.csv"), row.names = TRUE)

# Custom color palette for stromal clusters
col_group <- c("#FED439FF","#709AE1FF","#8A9197FF","#D2AF81FF","#FD7446FF","#D5E4A2FF",
               "#197EC0FF","#46732EFF","#71D0F5FF","#F79D1EFF","#1A9993FF","#FD8CC1FF","#7876B1FF")

# Predefined sample order for x-axis
sample_order <- c("C120","C20","C70","C90","C110","C130","C30","C40","C80","C50","C60","C100",
                  "C92","C122","C132","C32","C42","C52","C62","C102","C12","C94","C134","C44",
                  "C54","C64","C104","N44","N64","N94","N114","N134")

# Generate proportion stacked bar plot
p1 <- ggplot(cell.prop, aes(x = sample, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity", position = "fill", width = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(rep(1,4), 'cm'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = col_group) +
  scale_x_discrete(limits = sample_order)

# Save proportion plot
ggsave(file.path(out_dir, "03stromal_per_sample_cell_proportion.pdf"),
       plot = p1, units="in", dpi=300, width=10, height=3, device="pdf")


##===================== 2. Stacked bar plot: absolute cell number of stromal subcluster per sample =====================
# Convert count table to long format dataframe
cell.num2 <- as.data.frame(tab_count)
colnames(cell.num2) <- c("cluster","patient","cellnum")

# Export raw cell count table to csv
write.csv(cell.num2, file.path(out_dir, "stromal_stromal.num.csv"), row.names = FALSE)
write.csv(as.data.frame.matrix(cell_proportion), file.path(out_dir, "stromal_stromal.num2.csv"), row.names = TRUE)

# Generate absolute cell number stacked bar plot
p2 <- ggplot(cell.num2, aes(x = patient, y = cellnum, fill = cluster)) +
  geom_bar(stat="identity", width = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = unit(rep(1,4), 'cm'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  scale_fill_manual(values = col_group) +
  scale_x_discrete(limits = sample_order)

# Save cell count plot
ggsave(file.path(out_dir, "03stromal_per_sample_cell_num.pdf"),
       plot = p2, units="in", dpi=300, width=10, height=3, device="pdf")