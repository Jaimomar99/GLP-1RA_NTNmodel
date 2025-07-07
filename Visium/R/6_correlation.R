library(Seurat)
library(tidyverse)
library(parallel)
library(optparse)
library(psych)

parser <- OptionParser()

parser <- add_option(parser, c("-m", "--manifest_path", action="store", type="character", help="path to sample manifest csv"))
parser <- add_option(parser, c("-g", "--grouping", action="store", type="character", help="name of group column in manifest (ex. diagnosis)"))
args <- parse_args(parser)

manifest <- read.csv(args$manifest_path)

print("Creating list of Seurat objects)")
list_seurat <- lapply(1:dim(manifest)[1], FUN = function(x) {
  alias_id <- manifest[x, "alias_id"]
  seurat_obj <- readRDS(paste0(dirname(args$manifest_path), "/processed/by_sample/", alias_id, ".rds"))
  return(seurat_obj)
})

common_genes <- Reduce(intersect, lapply(list_seurat, rownames)) %>% sort()
conds <- unique(unlist(lapply(list_seurat, function(x) x[[args$grouping]])))


##create folders
out_dir <- paste0(dirname(args$manifest_path), "/processed/correlations/")
invisible(lapply(c("correlations_by_sample", "correlations_by_cond", "correlation_all"), function(x) {
  dir.create(file.path(out_dir))
  dir.create(file.path(out_dir, x))
}))


#filter objects based only common genes
list_seurat <- lapply(list_seurat, function(x) subset(x, features=common_genes))
gc() #liberate space

print("starting correlation by sample")

#correlation by sample
i = 0
invisible(mclapply(1:length(list_seurat), function(x){
  obj <- list_seurat[[x]]
  i <- i + 1
  message("Correlation analyzed:",i, '/', length(list_seurat))
  
  alias_id <- unique(obj$alias_id)
  
  df <- t(as.matrix(obj@assays$SCT@data))
  
  tmp <- psych::corr.test(df, method="spearman", use="pairwise.complete.obs")
  
  cor_matrix <- tmp$r
  pvalue_matrix <- tmp$p #only p, not p.adj
  
  write.csv(cor_matrix, paste0(out_dir,"correlations_by_sample/", alias_id, "_corr.csv"), row.names = TRUE )
  write.csv(pvalue_matrix, paste0(out_dir,"correlations_by_sample/", alias_id, "_pvalue.csv"), row.names = TRUE)
},mc.cores=16))

gc()


#correlation by cond
print("starting correlation by cond")
invisible(mclapply(conds, function(cond){
  i <- which(cond == conds)
  message("Conds analyzed:",i, '/', length(conds))
  
  cond_list <- Filter(function(y) cond %in% levels(pull(y[[args$grouping]])), list_seurat)
  
  expr_list <-  lapply(cond_list, function(x) t(as.matrix(x@assays$SCT@data)))
  
  expr_matrix <- do.call(rbind, expr_list)
  
  tmp <- psych::corr.test(expr_matrix, method="spearman", use="pairwise.complete.obs")
  
  cor_matrix <- tmp$r
  pvalue_matrix <- tmp$p
  write.csv(cor_matrix, paste0(out_dir,"correlations_by_cond/", cond, "_corr.csv" ) ,row.names = TRUE )
  write.csv(pvalue_matrix, paste0(out_dir,"correlations_by_cond/", cond, "_pvalue.csv"),row.names = TRUE)
},mc.cores=4))

gc()

#all_samples
print("starting correlation all samples")

expr_list <-  lapply(list_seurat, function(x) t(as.matrix(x@assays$SCT@data)))
expr_matrix <- do.call(rbind, expr_list)
tmp <- psych::corr.test(expr_matrix, method="spearman", use="pairwise.complete.obs")

cor_matrix <- tmp$r
pvalue_matrix <- tmp$p
write.csv(cor_matrix, paste0(out_dir,"correlation_all/", "corr.csv" ) ,row.names = TRUE )
write.csv(pvalue_matrix, paste0(out_dir,"correlation_all/", "pvalue.csv"),row.names = TRUE)

print("DONE")
