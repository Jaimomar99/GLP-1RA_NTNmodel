library(tidyverse)
library(Seurat)
library(optparse)

parser <- OptionParser()

parser <- add_option(parser, c("-m", "--manifest_path", action="store", type="character", help="path to sample manifest csv"))
parser <- add_option(parser, c("-o", "--out_dir", action="store", type="character", help="output subdirectory in processed"))
parser <- add_option(parser, c("-p", "--pipeline_path", action="store", type="character", help="path to pipeline to use"))
args <- parse_args(parser)
print(args)

print("Reading manifest ...")
source(paste0(args$pipeline_path, "/R/utils.R"))
manifest <- read.csv(args$manifest_path)
manifest$processed <- paste(dirname(args$manifest_path),"/processed/spaceranger", manifest$alias_id, "outs/", sep= "/")

n_samples <- nrow(manifest)

# Create output directories that might not exist
by_sample_out_dir <- paste0(dirname(args$manifest_path),"/processed/", args$out_dir)
report_out_dir <- paste0(dirname(args$manifest_path),"/reports/", args$out_dir)

if (!dir.exists(by_sample_out_dir)) {
  dir.create(by_sample_out_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(report_out_dir)) {
  dir.create(report_out_dir, recursive = TRUE, showWarnings = FALSE)
}

print("Running by sample normalization and clustering:")

list_seurat_summary <- lapply(1:n_samples, FUN = function(x) {
  current_sample <- manifest[x, ]
  print(current_sample)
  sample_alias <- current_sample$alias_id
  processed_path <- current_sample$processed
  filename <- list.files(processed_path, pattern = "filtered.*\\.h5$")
  
  print("Initiating seurat object")
  # create Seurat object for each sample
  seurat_obj <- Load10X_Spatial(data.dir = processed_path, filename = filename, slice = paste0("sample", sample_alias), use.names = TRUE)
  
  # add sample metadata to each object
  print("Adding sample metadata to seurat object")
  all_columns <- colnames(current_sample)
  all_columns = all_columns[!(all_columns == "fw_sample_id")]  # omit the fw sample id from metadata

  for (n in all_columns) {
    metadata_factor <- as.factor(current_sample[, n, drop = TRUE])
    seurat_obj <- AddMetaData(seurat_obj, metadata = metadata_factor, col.name = n) 
  }

  # summary stats b/f filter
  nspots <- ncol(seurat_obj)
  spot_median_umi <- median(seurat_obj[[]][, "nCount_Spatial"])
  spot_median_gene <- median(seurat_obj[[]][, "nFeature_Spatial"])

  print("Running QC") 
  seurat_obj <- filter_spatial(seurat_obj) # utils

  # summary stats after filter
  nspots_filter <- ncol(seurat_obj)

  spot_umi <- colSums(GetAssayData(seurat_obj, assay = "Spatial", slot = "counts"))  # re-calcualte post filter b/c Seurat doesn't update them
  spot_gene <- colSums(GetAssayData(seurat_obj, assay = "Spatial", slot = "counts") > 0)  # re-calcualte post filter b/c Seurat doesn't update them
  spot_median_umi_filter <- median(spot_umi)
  spot_median_gene_filter <- median(spot_gene)
  median_mito_filter <- median(seurat_obj[[]][, "percent_mito"])

  summary_stats_df <- data.frame(
    "sample_alias" = sample_alias,
    "total_spots" = nspots,
    "mediam_umi_per_spot" = spot_median_umi,
    "median_genes_per_spot" = spot_median_gene,
    "total_spots_filtered" = nspots_filter,
    "mediam_umi_per_spot_filtered" = spot_median_umi_filter,
    "median_genes_per_spot_filtered" = spot_median_gene_filter,
    "median_percent_mito_filtered" = median_mito_filter
  )

  print("Running normalization and clustering") 
  seurat_obj <- Normalization_clustering(seurat_obj, reduction="pca", sample=TRUE) #utils

  print("Writing seurat object to RDS")
  
  saveRDS(seurat_obj, file = paste0(by_sample_out_dir, "/", sample_alias, ".rds"))
  return(list("seurat_obj" = seurat_obj, "summary_stats" = summary_stats_df))
})

list_ST <- lapply(list_seurat_summary, function(x) x$seurat_obj)
list_summary_stats <- lapply(list_seurat_summary, function(x) x$summary_stats)

rm(list_seurat_summary)
gc()

gene_list <- lapply(1:n_samples, FUN = function(x) {
  seurat_obj <- list_ST[[x]]
  return(rownames(seurat_obj))
})

print("calculating gene freq")
# calculate gene detection frequency
keep.genes.combined <- as.data.frame(table(unlist(gene_list))) %>% dplyr::filter(Freq >= round(length(gene_list) * 0.2))
colnames(keep.genes.combined) <- c("gene_name", "freq")
keep.genes.combined$gene_name <- as.character(keep.genes.combined$gene_name)
write.csv(keep.genes.combined, file = paste0(report_out_dir, "/gene_detection_freq.csv"))

names(list_ST) <- paste0("sample", manifest$alias_id)

# print(list_ST)

print("Making raw tissue plot")
raw_plot <- plot_raw_data(list_ST)
ggsave(filename = paste0(report_out_dir, "/raw_tissue.png"), plot = raw_plot)

# merge summary stats list
summary_stats_combined <- do.call(rbind, list_summary_stats)
rownames(summary_stats_combined) <- manifest$alias_id
write.csv(summary_stats_combined, file = paste0(report_out_dir, "/sample_summary_stats.csv"), row.names = FALSE)
