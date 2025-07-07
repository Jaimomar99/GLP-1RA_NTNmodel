library(optparse)
library(fs)

parser <- OptionParser()
parser <- add_option(parser, c("-n", "--dataset_name", action="store", type="character", help="dataset name", default=NULL))
parser <- add_option(parser, c("-d", "--dir_path", action="store", type="character", help="absolute path to directory where
    dataset should be set up", default=NULL))

args <- parse_args(parser)

# Check if project name is provided
if (is.null(args$dataset_name)) {
    stop("Please provide a name for the project directory using the --dataset_name option.")
}

if (is.null(args$dir_path)) {
    stop("Please provide a path to directory where project should be created")
}

# Create project directory if it doesn't already exist
project_dir <- path(args$dir_path, args$dataset_name)
dir_create(project_dir)

# Set working directory to project dir
setwd(project_dir)

# Create project subdirectories if they don't exist
dir_create("raw", "fastq")
dir_create("raw", "images")
dir_create("raw", "manual_align_image")
dir_create("processed", "spaceranger")
dir_create("processed", "by_sample")
dir_create("processed", "integrated")
dir_create("processed", "deconvolution")
dir_create("processed", "loupe_browser")
dir_create("processed", "differential_expression")
dir_create("processed", "loupe_browser")
dir_create("processed", "ultimate_processed", "by_sample")
dir_create("processed", "ultimate_processed", "integrated")
dir_create("reports")
dir_create("notebooks")
dir_create("workflow")
dir_create("workflow", "logs")
