suppressWarnings(library(Seurat))
suppressWarnings(library(optparse))
suppressWarnings(library(tidyverse))

parser <- OptionParser()

parser <- add_option(parser, c("-m", "--manifest_path", action="store", type="character", help="path to sample manifest csv"))
parser <- add_option(parser, c("-i", "--integrated_out", action="store", type="character", help="path to results integrated object"))
parser <- add_option(parser, c("-o", "--out_dir", action="store", type="character", help="output subdirectory in processed"))
parser <- add_option(parser, c("-p", "--pipeline_path", action="store", type="character", help="path to pipeline to use"))

args <- parse_args(parser)

print(args)
source(paste0(args$pipeline_path, "/R/utils.R"))

integrated <- readRDS(paste0(dirname(args$manifest_path),"/processed/", args$integrated_out , "/harmonized.rds"))

#calculate DE genes per cluster in the integrated object
de_markers <- DE_function(integrated, "clustering_harmony", assay="SCT")

# create out directory if it doesn't exist
out_dir <- paste0(dirname(args$manifest_path),"/processed/differential_expression/", args$integrated_out, "/")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

write.csv(de_markers, paste0(out_dir, "/de_markers.csv"), row.names=FALSE)

print("DE markers calculated and saved")