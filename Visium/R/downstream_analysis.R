suppressWarnings(library(Seurat))
suppressWarnings(library(optparse))
suppressWarnings(library(tidyverse))

parser <- OptionParser()

parser <- add_option(parser, c("-m", "--manifest_path", action="store", type="character", help="path to sample manifest csv"))
parser <- add_option(parser, c("-o", "--deconvolute_out", action="store", type="character", help="path to results deconvolution"))
parser <- add_option(parser, c("-s", "--deconvolute_nsplits", action="store", type="int", help="number of nsplits in deconvolution"))
parser <- add_option(parser, c("-i", "--integrated_out", action="store", type="character", help="path to results integrated object"))
parser <- add_option(parser, c("-b", "--by_sample_out", action="store", type="character", help="path to results by sample object"))
parser <- add_option(parser, c("-u", "--ultimate_processed_out", action="store", type="character", help="path to ultimate folder results"))
parser <- add_option(parser, c("-l", "--loupe_browser_out", action="store", type="character", help="path to loupe browser"))
parser <- add_option(parser, c("-p", "--pipeline_path", action="store", type="character", help="path to pipeline to use"))
parser <- add_option(parser, c("-f", "--feature_df_path", action="store", type="character", help="path to csv file with feature info; gene_id column with ensembl ids and gene_symbol column"))

args <- parse_args(parser)

print(args)
source(paste0(args$pipeline_path, "/R/utils.R"))


print("Reading manifest ...")
manifest <- read.csv(args$manifest_path)
abund_path <- paste(dirname(args$manifest_path), "/processed/", args$deconvolute_out, sep= "/")

# check features dataframe; set to NULL if length = 0
if (length(args$feature_df_path) == 0) {
  args$feature_df_path = NULL
}

#read abundance
if (args$deconvolute_nsplits > 1){
    pattern <- paste0("^split_(", paste0(0:args$deconvolute_nsplits, collapse = "|"), ")")

    files <- list.files(abund_path, full.name=TRUE, pattern = pattern)
    
    abundances <- lapply(files, function(x){
        read.csv(paste(x, "cell_abundance_q5.csv", sep = "/"), row.names=1)
    })
    
    # Combine all data frames into one
    abund_df <- do.call(rbind, abundances)
    print("deconvolution abundance dataframes merged")
} else {
    abund_df <- read.csv(paste(abund_path, "cell_abundance_q5.csv", sep = "/"), row.names=1)
    print("deconvolution abundance read")
}

#read integrated object
integrated_path <- paste0(dirname(args$manifest_path),"/processed/", args$integrated_out , "/harmonized.rds")
integrated <- readRDS(integrated_path)

abund_df <- abund_df[match(rownames(integrated@meta.data), rownames(abund_df)), ]  # reorder abundance rows to match integrated metadata
integrated@meta.data <- cbind(integrated@meta.data, abund_df)
print("Integrated read and merged with deconvolution proportions")

#clustering by abundance
integrated<- infer_CC(integrated, abund_df)
print("Compositional clusters calculated")

# run log1p counts per 10K normalization on combined/integrated counts
print("Calculating log1p counts per 10K")
integrated <- NormalizeData(integrated, assay = "Spatial", normalization.method = "LogNormalize", scale.factor = 10000)  # using defaults
integrated <- ScaleData(integrated, assay = "Spatial")  # using defaults

#create directories
lapply(c("rds/", "meta/", "counts/", "SCT_counts/", "logCP10K/", "tissue_positions/"), function(d) {
  dir.create(paste0(dirname(args$manifest_path), "/processed/", args$ultimate_processed_out, "/by_sample/", d), recursive = TRUE)
})

#loop by sample
n_samples <- dim(manifest)[1]
cell_types <- abund_df %>% colnames()

lapply(1:n_samples, FUN = function(x) {
    current_sample <- manifest[x, ]
    sample_alias <- current_sample$alias_id
    
    seurat_obj_sample <- readRDS(paste0(dirname(args$manifest_path), "/processed/", args$by_sample_out, "/", sample_alias, ".rds"))
    # run log1p counts per 10K normalization on combined/integrated counts
    seurat_obj_sample <- NormalizeData(seurat_obj_sample, assay = "Spatial", normalization.method = "LogNormalize", scale.factor = 10000)  # using defaults
    seurat_obj_sample <- ScaleData(seurat_obj_sample, assay = "Spatial")  # using defaults
    
    #adding to each sample, clustering by integration, clustering by composition and deconvolution 
    subset_integrated <- integrated@meta.data[, c("alias_id", "clustering_harmony", colnames(abund_df), "Compositional_Cluster")] %>% 
      filter(alias_id == sample_alias) %>% 
      mutate(seq = rownames(.)) %>% 
      mutate(seq = gsub("_\\d+$", "", seq))
    
    rownames(subset_integrated) <- subset_integrated$seq
    seurat_obj_sample@meta.data <- merge(seurat_obj_sample@meta.data, subset_integrated, by = 0, all.x = TRUE, suffixes = "") %>% 
      column_to_rownames("Row.names") %>% 
      select(!c("seq", "alias_idNA"))
    
    #save rds, counts, norm & meta
    output_path <- paste0(dirname(args$manifest_path), "/processed/", args$ultimate_processed_out, "/by_sample/")
    write_seurat_spatial_coord(seurat_obj_sample, paste0(output_path, "/tissue_positions/"), integrated = FALSE)
    save_counts_norm_meta(seurat_obj_sample, sample_alias, output_path, feature_df_path = args$feature_df_path)
    
    # save dimensionality reduction
    lapply(c("pca", "umap"), function(r) {
      write_dim_reduction(seurat_obj_sample, paste0(output_path, "/dim_reduction/", sample_alias, "/"), dim_reduction = r)
    })
    print(paste("sample: ", sample_alias, "counts, meta, norm & rds saved", sep = " "))
    

    #create loupe browser files
    invisible(lapply(c("Compositional_Cluster", "clustering_by_sample", "clustering_harmony"), function(column){
      column_path <- paste0(dirname(args$manifest_path), "/processed/", args$loupe_browser_out, "/", column)
      if (!dir.exists(column_path)) {dir.create(column_path, recursive = TRUE)}    # If the folder does not exist, create it
      
      #save loupe files
      seurat_obj_sample@meta.data %>% 
        dplyr::select(column) %>% 
        rownames_to_column("Barcodes") %>% 
        write_csv(paste0(column_path, "/", sample_alias, ".csv"))
    }))
    print(paste0("saving sample:", sample_alias, " loupe browsers"))

    #save cell types proportions plots
    reports_path <- paste0(dirname(args$manifest_path), "/reports/", args$deconvolute_out, "/", sample_alias, "/" )
    if (!dir.exists(reports_path)) {dir.create(reports_path, recursive = TRUE)}
    
    lapply(cell_types, function(x){
      p <- SpatialFeaturePlot(seurat_obj_sample,
                features = x,
                 image.alpha = 0.5,
                 crop=FALSE,
                 pt.size.factor = 1.2             
                )
  
      ggsave(paste0(reports_path, "/", x, ".png"), plot = p, width = 18, height = 15)
    })
    
    # save piechart proportion
    p <- plotSpatialScatterpie(seurat_obj_sample, abund_df, pie_scale=0.4)
    ggsave(paste0(reports_path, "/pie_chart.png"), plot = p, width = 18, height = 15)
    print(paste0("cell types plots of sample:", sample_alias, " saved"))
})

#saving integrated object
lapply(c("rds/", "meta/", "counts/", "SCT_counts/", "logCP10K/", "tissue_positions/", "dim_reduction"), function(d) {
  dir.create(paste0(dirname(args$manifest_path), "/processed/", args$ultimate_processed_out, "/integrated/", d), recursive = TRUE)
})

output_path <- paste0(dirname(args$manifest_path), "/processed/", args$ultimate_processed_out, "/integrated")
write_seurat_spatial_coord(integrated, paste0(output_path, "/tissue_positions/"), integrated = TRUE)
save_counts_norm_meta(integrated, "integrated", output_path, feature_df_path = args$feature_df_path)
lapply(c("pca", "umap", "harmony"), function(r) {
  write_dim_reduction(integrated, paste0(output_path, "/dim_reduction/"), dim_reduction = r)
})

p <- Composition_cluster_enrichedCelltypes(integrated, abund_df)

ggsave(paste0(dirname(args$manifest_path), "/reports/", args$deconvolute_out, "/Deconvolution_proportions.png"), plot = p, width = 18, height = 15)

print("integrated object saved script finished")