## create clusters and groups at different resolutions and using different PCs in the first place

##--------------------------------------------------------------------------------------------------##
rm(list = ls())
set.seed(1234)
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
library(Seurat, lib.loc = "/ceph/package/c7/R-cbrg/current/4.2.0")
library(optparse, lib.loc = "/ceph/package/c7/R-cbrg/current/4.2.0")
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
option_list = list(
  make_option(c("--nPCs"), type="character", default=NA,
              help="Number of principal components to use for clustering", metavar="numeric"))
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
##--------------------------------------------------------------------------------------------------##
  
##--------------------------------------------------------------------------------------------------##
robject_folder <- "~/projects/lighting/data/robjects/" ## CCB cluster
genome <- 'dmel649ChrimsonV2'
input_objs <- paste0('KB-', genome)
robject_folder <- paste0(robject_folder, input_objs, '/')
robject_folder2 <- paste0(robject_folder, 'clustering_nPCs/')
system(paste0("mkdir -p ", robject_folder))
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## functions
source("../helper_functions.R")
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Analyse the integrated dataset
combined_dataset <- readRDS(paste0(robject_folder, 'combined_dataset.RDS'))
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## re-calculate groups for integrate dataset -- I don't know why there are clusters already, possibly from the individual datasets?
DefaultAssay(combined_dataset) <- 'integrated'

## remove any metadata from previously-computed clusters
try({
combined_dataset[[]][,grepl('integrated_snn_res', colnames(combined_dataset[[]]))] <- NULL
})

## compute
combined_dataset <- FindNeighbors(combined_dataset, dims = 1:opt$nPCs)
combined_dataset <- FindClusters(combined_dataset, dims = 1:opt$nPCs)

res_vect <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5)

for(res_it in res_vect){
  combined_dataset <- FindClusters(combined_dataset, dims = 1:opt$nPCs, resolution = res_it)
}

## save only metadata with the clusters
saveRDS(combined_dataset[[]][,grepl('integrated_snn_res', colnames(combined_dataset[[]]))],
        paste0(robject_folder2, 'clusters_nPCs', opt$nPCs, '.RDS'))
##--------------------------------------------------------------------------------------------------##
