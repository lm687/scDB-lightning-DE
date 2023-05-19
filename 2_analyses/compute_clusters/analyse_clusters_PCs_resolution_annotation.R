rm(list = ls())

library(fossil)
library(mcclust)
library(ggpointdensity)
library(scico)
library(ComplexHeatmap)

##-----------------------------------------------------------------------------------------------------------------##
robject_folder <- "~/projects/lighting/data/robjects/" ## CCB cluster
genome <- 'dmel649ChrimsonV2'
input_objs <- paste0('KB-', genome)
robject_folder <- paste0(robject_folder, input_objs, '/')
folder_results <- "~/projects/lighting/3_results/" ## CCB cluster

input_fles_clusters <- list.files("/home/l/lmorrill/projects/lighting/data/robjects/KB-dmel649ChrimsonV2/clustering_nPCs/", full.names = T)

files_clusters <- lapply(input_fles_clusters, readRDS)

##-----------------------------------------------------------------------------------------------------------------##

## get expression matrix
combined_dataset <- readRDS(paste0(robject_folder, 'combined_dataset.RDS'))

files_clusters

## findmarkers should be done on the RNA, not integrated, dataset
DefaultAssay(combined_dataset) <- 'RNA'
combined_dataset@assays$RNA@data

## for each of the clustering combinations, find the most common markers
combined_dataset
sapply(1:length(files_clusters), function(npcs){
  sapply(1:ncol(files_clusters[[1]]), function(res){
    cat('npcs =', npcs, '; res =', res, '\n')
    combined_dataset_mod <- combined_dataset
    # combined_dataset_mod$seurat_clusters <- files_clusters[[npcs]][,res]
    Idents(combined_dataset_mod) <- files_clusters[[npcs]][,res]
    mrkrs <- Seurat::FindAllMarkers(combined_dataset_mod)
    saveRDS(mrkrs, paste0("/home/l/lmorrill/projects/lighting/data/robjects/KB-dmel649ChrimsonV2/clustering_nPCs_markers/markers_res", res, "_npcs", npcs, ".RDS"))
  })
})

