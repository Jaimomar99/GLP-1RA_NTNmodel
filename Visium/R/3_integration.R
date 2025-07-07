library(tidyverse)
library(Seurat)
library(harmony)
library(optparse)
library(SeuratDisk)
library(stringr)

parser <- OptionParser()

parser <- add_option(parser, c("-m", "--manifest_path", action="store", type="character", help="path to sample manifest csv"))
parser <- add_option(parser, c("-s", "--by_sample_subdir", action="store", type="character", help="by sample sub directory"))
parser <- add_option(parser, c("-o", "--out_dir", action="store", type="character", help="output subdirectory in processed"))
parser <- add_option(parser, c("-g", "--grouping", action="store", type="character", help="name of group column in manifest (ex. diagnosis)"))
parser <- add_option(parser, c("-b", "--batch_correction_var", action="store", type="character", help="list of variables that harmony should correct for"))
parser <- add_option(parser, c("-p", "--pipeline_path", action="store", type="character", help="path to pipeline to use"))
args <- parse_args(parser)

print(args)
source(paste0(args$pipeline_path, "/R/utils.R"))

print("Reading manifest ...")
manifest <- read.csv(args$manifest_path)

manifest$by_sample <- paste0(dirname(args$manifest_path),"/processed/", args$by_sample_subdir, "/", manifest$alias_id, ".rds")
manifest$processed <- paste(dirname(args$manifest_path),"/processed/spaceranger", manifest$alias_id, "outs/", sep= "/")
n_samples <- nrow(manifest)

# create output directories that don't exist
integration_out_dir <- paste0(dirname(args$manifest_path),"/processed/", args$out_dir)
report_out_dir <- paste0(dirname(args$manifest_path),"/reports/", args$out_dir)

if (!dir.exists(integration_out_dir)) {
  dir.create(integration_out_dir, recursive = TRUE, showWarnings = FALSE)
}

if (!dir.exists(report_out_dir)) {
  dir.create(report_out_dir, recursive = TRUE, showWarnings = FALSE)
}

batch_correction_var = unlist(str_split(args$batch_correction_var, "\\s*,\\s*"))
print(batch_correction_var)

# assuming that gene detection csv exists
keep.genes.combined <- read.csv(paste0(dirname(args$manifest_path), "/reports/", args$by_sample_subdir, "/gene_detection_freq.csv"), stringsAsFactor = FALSE)
keep.genes <- keep.genes.combined$gene_name

print("Creating list of Seurat objects)")
list_ST <- lapply(1:n_samples, FUN = function(x) {
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

  # only filter spots
  seurat_obj <- filter_spatial(seurat_obj, min_spot_fraction=0.0)
  seurat_obj <- seurat_obj[rownames(seurat_obj) %in% keep.genes, ]

  seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = TRUE, vst.flavor = "vst2")
  DefaultAssay(seurat_obj) <- "SCT"

  return(seurat_obj)
})

print("Running harmony correction")
harmonized_obj <- Harmony_integration(list_ST, batch_cond_vector = batch_correction_var)

harmonized_obj <- Normalization_clustering(harmonized_obj, reduction="harmony", sample=FALSE)

print("Making plots")
QC <- QC_plot(
    harmonized_ST = harmonized_obj,
    subject_id_column = "alias_id",
    group_column = args$grouping,
    count_column = "nCount_Spatial",
    feature_column = "nFeature_Spatial",
    percent.mt = "percent_mito",
    spots = NULL
)

ggsave(filename = paste0(report_out_dir, "/QC.png"), plot = QC)

harmony_plot <- harmony_checkpoint_plot(harmonized_obj, batch_cond_vector = batch_correction_var)
ggsave(filename = paste0(report_out_dir, "/harmony_plot.png"), plot = harmony_plot)

print("Running marker prep")
harmonized_obj <- PrepSCTFindMarkers(harmonized_obj)

print("Writing harmony corrected seurat object to rds")
#saving objects
saveRDS(harmonized_obj, file = paste0(integration_out_dir, "/harmonized.rds"))

total_spots <- ncol(harmonized_obj)
spot_umi <- colSums(GetAssayData(harmonized_obj, assay = "Spatial", slot = "counts"))  # re-calcualte post filter b/c Seurat doesn't update them
spot_gene <- colSums(GetAssayData(harmonized_obj, assay = "Spatial", slot = "counts") > 0)  # re-calcualte post filter b/c Seurat doesn't update them
spot_median_umi <- median(spot_umi)
spot_median_gene <- median(spot_gene)
median_percent_mito <- median(harmonized_obj$percent_mito)

summary_stats_df <- data.frame(
  "total_spots" = total_spots,
  "mediam_umi_per_spot" = spot_median_umi,
  "median_genes_per_spot" = spot_median_gene,
  "median_percent_mito" = median_percent_mito
) 

write.csv(summary_stats_df, file = paste0(report_out_dir, "/dataset_summary_stats.csv"))

# getting errors with Diet Seurat when trying to remove data
# diet_obj <- DietSeurat(
#   harmonized_obj,
#   counts=TRUE,
#   data=FALSE,
#   assays = "SCT",
#   dimreducs = FALSE,
#   misc=FALSE)

# manually created diet_obj
diet_obj <- CreateSeuratObject(counts = harmonized_obj@assays$Spatial@counts, meta.data = harmonized_obj@meta.data)

print("writing harmony corrected suerat object to h5Seurat")
SaveH5Seurat(diet_obj, file = paste0(integration_out_dir, "/combined.h5Seurat"))

print("converting seurat to anndata h5ad")
Convert(paste0(integration_out_dir, "/combined.h5Seurat"), dest = "h5ad",  overwrite = TRUE)
file.remove(paste0(integration_out_dir, "/combined.h5Seurat"))

# TO DO: move this to it's own script
# print("Running cluster differential expression")
# de_markers <- DE_function(harmonized_obj, "clustering_harmony", assay="SCT")
# 
# write_csv(de_markers, paste0(integration_out_dir,"differential_expression/de_markers_clusters.csv"))
