library(infercnv)
library(tidyverse)

#==================== Parameter Section [Modify paths here] ====================
# Directory for kmeans_df csv files
kmeans_dir <- "./test/infercnv/00kmeans_df"
# Root directory of infercnv outputs
infercnv_root <- "./test/infercnv/result/epi/"
# Output directory for expression csv files
out_expr_dir <- "./test/infercnv/02expr_z"
# Total sample number
sample_num <- 9
# Whether perform z-score normalization on inferCNV matrix
do_zscore <- TRUE
# Whether enable gene filtering
do_filter_gene <- FALSE
#=================================================================

#==== 1. Read and merge all kmeans_df files ====
setwd(kmeans_dir)
# Search csv files containing "kmeans_df"
kmeans_files <- list.files(pattern = "kmeans_df.*csv$", full.names = FALSE)
if(length(kmeans_files) == 0) stop("ERROR: No csv files starting with kmeans_df found in kmeans_dir!")
# Only keep first N samples
kmeans_files <- head(kmeans_files, sample_num)

# Batch read files, use bind_rows instead of rbind loop for higher efficiency
newdatat <- map_dfr(kmeans_files, function(f){
  df <- read.csv(f)
  return(df)
})
# Rename columns
colnames(newdatat) <- c("cell_name", "kmeans_class", "group")
cat("Total annotated cells after merging:", nrow(newdatat), "\n")
head(newdatat)
# Save combined annotation table as backup
saveRDS(newdatat, file.path(kmeans_dir, "combined_kmeans_anno.rds"))

#==== 2. Batch load all infercnv objects ====
# Find folders starting with "C"
folders <- list.files(infercnv_root, pattern = "^C", full.names = FALSE)
folders <- head(folders, sample_num)
cat("Samples to process:", paste(folders, collapse = ","), "\n")

# Generate full paths to run.final.infercnv_obj
folders_file <- map_chr(folders, function(sample_name){
  file.path(infercnv_root, sample_name, "run.final.infercnv_obj")
})
# Check existence of each infercnv file
file_exist <- file.exists(folders_file)
if(!all(file_exist)){
  cat("WARNING: Some infercnv files missing:", folders_file[!file_exist], "\n")
  folders_file <- folders_file[file_exist]
  folders <- folders[file_exist]
}

# Function to read infercnv object
# path: path to run.final.infercnv_obj
# filter: whether filter genes
# zscore: whether perform row-wise z-score normalization on genes
read_infercnv_obj <- function(path, filter=FALSE, zscore=TRUE) {
  
  if(!file.exists(path)){
    warning(paste("File not found:", path))
    return(NULL)
  }
  
  # Extract sample ID from folder name
  sample_folder <- basename(dirname(path))
  sample_split <- str_split(sample_folder, "_")[[1]]
  sample_id <- tail(sample_split, 1)
  
  infercnv_obj <- readRDS(path)
  gene_order <- as.data.frame(infercnv_obj@gene_order)
  
  expr <- infercnv_obj@expr.data # matrix: genes × cells
  
  if(zscore) {
    # Row-wise z-score for each gene, avoid division by zero when sd=0
    expr <- t(apply(expr, 1, function(y){
      mu <- mean(y)
      s <- sd(y)
      if(s == 0){
        return(rep(0, length(y)))
      }else{
        return((y - mu)/s)
      }
    }))
  }
  
  if (filter) {
    # Filter genes: number of values >3 exceeds 1/10 of total cells
    expr <- expr[apply(expr, 1, function(row) {
      sum(row > 3) > (length(row) / 10)
    }), ]
  }
  
  # Retain genes present in both expr matrix and gene_order
  keep_genes <- intersect(rownames(expr), rownames(gene_order))
  expr_filtered <- expr[keep_genes, ]
  
  cat(sprintf("Sample %s, n_genes:%d, n_cells:%d\n", sample_id, nrow(expr_filtered), ncol(expr_filtered)))
  return(expr_filtered)
}

# Read infercnv objects with tryCatch for error handling
sceList <- map(folders_file, function(fp){
  tryCatch({
    read_infercnv_obj(fp, filter = do_filter_gene, zscore = do_zscore)
  }, error = function(e){
    cat("Failed to read:", fp, "\n", e$message, "\n")
    return(NULL)
  })
})
names(sceList) <- folders
# Remove NULL entries caused by read failure
sceList <- compact(sceList)
folders_valid <- names(sceList)

#==== 3. Subset cells matched with kmeans annotation and export csv ====
setwd(out_expr_dir)
for (i in seq_along(sceList)) {
  sample_name <- folders_valid[i]
  mat <- sceList[[i]]
  # Find common cell names
  common_cells <- intersect(colnames(mat), newdatat$cell_name)
  cat(sprintf("Sample %s, total cells in infercnv:%d, matched annotated cells:%d\n",
              sample_name, ncol(mat), length(common_cells)))
  
  if(length(common_cells) == 0){
    warning(paste0(sample_name, " has zero matched cells, skip this sample"))
    next
  }
  mat_sub <- mat[, common_cells, drop = FALSE]
  
  out_csv <- paste0(sample_name, ".csv")
  write.csv(mat_sub, file = out_csv, row.names = TRUE)
  cat("Output file:", out_csv, "\n")
  
  # Free memory after each iteration
  gc()
}
cat("All samples export finished!\n")