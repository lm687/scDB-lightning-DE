# ## analyse light dataset with bioconductor packages

#--------------------------------------------------------------------------------#

rm(list = ls())
set.seed(1234)

library(Seurat)
library(ggplot2)
library(cowplot)
library(plotly)
library(reshape2)
library(gridExtra)
library(DropletUtils)
library(tidyr)
library(dplyr)
theme_set(theme_cowplot())

#--------------------------------------------------------------------------------#

local <- T
if(local){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  # folder_input <- '../../github-repo-lightning-DE/3_results_local/dmel649Chrimson/' ## not all files are available
  folder_input <- '../../github-repo-lightning-DE/3_results_local/dmel649CORRECTED//' ## not all files are available
}else{
  ## CCB cluster
  setwd("/t1-data/project/cncb/shared/proj002/analyses/2020/10kCells/")
  folder_input <- './' 
}

source("helper_functions.R")

#--------------------------------------------------------------------------------#

samples <- apply(expand.grid(paste0('G', 1:4), paste0('_rep', 1:2)), 1, paste0, collapse='')
samples

input_files <- lapply(samples, function(rep_it){
  cat(rep_it, '\n')
  # filtered_folder <- paste0(folder_input, rep_it, '/', rep_it, '/outs/filtered_feature_bc_matrix')
  raw_folder <- paste0(folder_input, rep_it, '/', rep_it, '/outs/raw_feature_bc_matrix')
  # print(raw_folder)
  DropletUtils::read10xCounts(raw_folder)
  # DropletUtils::read10xCounts(filtered_folder)
})
names(input_files) <- samples

input_files_single_obj <- DropletUtils::read10xCounts(paste0(folder_input, samples, '/', samples, '/outs/raw_feature_bc_matrix'))

#--------------------------------------------------------------------------------#

# input_files_Seurat <-  lapply(1:length(samples), function(i){
#   CreateSeuratObject(input_files[[i]], project = samples[i],
#                                     min.cells = 3, min.features = 200)
# })
# names(input_files_Seurat) <- samples

# input_files_Seurat <-  lapply(1:length(samples), function(i){
#   input_files_Seurat[[i]]$percent.mt <- PercentageFeatureSet(input_files_Seurat[[i]], pattern = "^mt:")
# })

#--------------------------------------------------------------------------------#

## change gene names
name_conversion_file <- give_name_conversion_file(path_tsv_genes = '/Users/lenamorrill/Documents/projects/general/genes/fb_synonym_fb_2022_06_cut.tsv')

for(i in 1:length(input_files)){
  rownames(input_files[[i]]) <- convert_FB_to_name(rownames(input_files[[i]]))
}

input_files$G1_rep1

rownames(input_files_single_obj) <- convert_FB_to_name(rownames(input_files_single_obj))

#--------------------------------------------------------------------------------#

# # COPYFEB2019_G1_rep1 <- FEB2019_G1_rep1
# # rep_it= samples[1]
# # for(rep_it in samples){
# #   assign(x = paste0('COPYFEB2019_', rep_it), pos = 'percent.mt',
# #          PercentageFeatureSet(get(paste0('FEB2019_', rep_it)), pattern = "^mt:"))
# #   # assign(x = paste0('COPYFEB2019_', rep_it, "[['percent.mt']]"), 
# #   #        PercentageFeatureSet(get(paste0('FEB2019_', rep_it)), pattern = "^mt:"))
# #   assign(x = paste0('COPYFEB2019_', rep_it, "@meta.data"), 
# #          cbind.data.frame(get(paste0('COPYFEB2019_', rep_it))@meta.data,
# #                           percent.mt=PercentageFeatureSet(get(paste0('FEB2019_', rep_it)), pattern = "^mt:")[,1])) ## not enough
# #   get(paste0('COPYFEB2019_', rep_it))[['percent.mt']] <- PercentageFeatureSet(get(paste0('FEB2019_', rep_it)), pattern = "^mt:") ## doesn't work
# #   # assign(x = paste0('FEB2019_', rep_it), 
# #   #        list(paste0('FEB2019_', rep_it, collapse=''), percent.mt=PercentageFeatureSet(get(paste0('FEB2019_', rep_it)), pattern = "^mt:")))
# # }
# 
# FEB2019_G1_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G1_rep1, pattern = "^mt:")
# FEB2019_G1_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G1_rep2, pattern = "^mt:")
# FEB2019_G2_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G2_rep1, pattern = "^mt:")
# FEB2019_G2_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G2_rep2, pattern = "^mt:")
# FEB2019_G3_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G3_rep1, pattern = "^mt:")
# FEB2019_G3_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G3_rep2, pattern = "^mt:")
# FEB2019_G4_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G4_rep1, pattern = "^mt:")
# FEB2019_G4_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G4_rep2, pattern = "^mt:")
# 
# do.call('grid.arrange', lapply(samples, function(rep_it) FeatureScatter(get(paste0('FEB2019_', rep_it)),
#                                                                         feature1 = "nCount_RNA", feature2 = "percent.mt")))
# 
# do.call('grid.arrange', lapply(samples, function(rep_it) FeatureScatter(get(paste0('FEB2019_', rep_it)),
#                                                                         feature1 = "nCount_RNA", feature2 = "nFeature_RNA")))

# barcoderanks_out <- lapply(input_files, function(rep_it) DropletUtils::barcodeRanks(m = rep_it))
# plot(as(barcoderanks_out$G1_rep1[,1:2], 'matrix'))

#--------------------------------------------------------------------------------#

# ## Empty droplets, doublets
emptydrops <- lapply(input_files, function(rep_it) DropletUtils::emptyDrops(m = assay(rep_it)))
emptydrops$G1_rep1

# swappeddrops <- lapply(input_files, function(rep_it) DropletUtils::swappedDrops(assay(rep_it)))

## QC
perCellQCMetrics_out <- lapply(input_files, function(i) scuttle::perCellQCMetrics(assay(i)))
perFeatureQCMetrics_out <- lapply(input_files, function(i) scuttle::perFeatureQCMetrics(assay(i)))

par(mfrow=c(2,4))
sapply(1:length(perCellQCMetrics_out), function(i){
  plot(log(as(perCellQCMetrics_out[[i]], 'matrix')),
       main=names(perCellQCMetrics_out)[i])
})

par(mfrow=c(2,4))
sapply(1:length(perFeatureQCMetrics_out), function(i){
  plot(log(as(perFeatureQCMetrics_out[[i]], 'matrix')),
       main=names(perFeatureQCMetrics_out)[i])
})


#--------------------------------------------------------------------------------#


## add percentage of mt, ncounts, nfeatures
is.mito <- lapply(input_files, function(i) grepl("^mt:", rownames(rowData(i))))
sapply(is.mito, mean)
for(i in 1:length(input_files)){
  input_files[[i]] <- scuttle::addPerCellQC(input_files[[i]], percent_top = 50,
                                            subsets = list(MT=is.mito[[i]]))
}


# plot1 <- lapply(input_files, function(i) colData(i) %>%
#   as_tibble() %>% 
#   ggplot() +
#   geom_violin(aes(Sample, sum)) +
#   labs(x = "Total UMI", y = "Value")+scale_y_continuous(trans = "log2"))
# plot2 <-  lapply(input_files, function(i) colData(i) %>%
#   as_tibble() %>% 
#   ggplot() +
#   geom_violin(aes(Sample, detected)) +
#   labs(x = "Genes detected", y = "Value")+scale_y_continuous(trans = "log2"))
# plot3 <- lapply(input_files, function(i) colData(i) %>%
#   as_tibble() %>% 
#   ggplot() +
#   geom_violin(aes(Sample, subsets_MT_percent)) +
#   labs(x = "Percentage mitochondrial", y = "Value")+scale_y_continuous(trans = "log2"))
# 
# cowplot::plot_grid(plot1)
# cowplot::plot_grid(plot2)
# do.call('grid.arrange', plot2)
# cowplot::plot_grid(plot3)
qc_violin <- (melt(lapply(input_files, function(i) data.frame(colData(i)[,c('sum', 'detected', 'subsets_MT_percent')]))))
ggplot(qc_violin, aes(x=L1, y=value))+geom_violin()+facet_wrap(.~L1)

plot(density(log(1+input_files[[1]]$'sum')))


#--------------------------------------------------------------------------------#

## Additional filtering

# ## Filter by other means
input_files <- lapply(input_files, function(i) i[, i$sum > 4500 & i$subsets_MT_percent < 15 & i$detected > 1500])

#--------------------------------------------------------------------------------#

## Normalisation
input_files <- lapply(input_files, scuttle::logNormCounts)

library(DelayedMatrixStats)

plot_MeanVar <- function(obj, count_name="counts"){
  x <- DelayedArray(assay(obj, count_name))
  plot_data <- tibble(
    mean = DelayedMatrixStats::rowMeans2(x),
    variance = DelayedMatrixStats::rowVars(x)
  )
  plot_data
  ggplot(plot_data, aes(mean, variance)) +
    geom_point()
}
do.call('grid.arrange', lapply(input_files, plot_MeanVar))
do.call('grid.arrange', lapply(input_files, plot_MeanVar, count_name = "logcounts"))

#--------------------------------------------------------------------------------#
## integration
input_files_single_obj <- scuttle::logNormCounts(input_files_single_obj)
input_files_single_obj

## see if there are batch effects


#--------------------------------------------------------------------------------#
## Feature selection
library(scran)
featuresec <- lapply(input_files, scran::modelGeneVar)
featuresec

do.call('grid.arrange', lapply(featuresec, function(i) ggplot(as_tibble(i)) +
  geom_point(aes(mean, total), color = "black") +
  geom_point(aes(mean, bio), color = "blue") +
  geom_point(aes(mean, tech), color = "red")))

TopHVGs <- lapply(featuresec, scran::getTopHVGs)

do.call('grid.arrange', lapply(1:length(featuresec), function(i) featuresec[[i]] %>%
  as_tibble() %>%
  mutate(
    gene_id = rownames(featuresec[[i]]),
    hvg = gene_id %in% TopHVGs[[i]]
  ) %>%
  ggplot() +
  geom_point(aes(mean, bio, color = hvg))))

#--------------------------------------------------------------------------------#
for(i in 1:length(input_files)){
  input_files[[i]] <- scater::runPCA(input_files[[i]], subset_row=TopHVGs[[i]])
}

for(i in 1:length(input_files)){
  input_files[[i]] <- scater::runUMAP(input_files[[i]], dimred = 'PCA', external_neighbors=TRUE)
}

for(i in 1:length(input_files)){
  input_files[[i]] <-  scater::runTSNE(input_files[[i]], dimred = 'PCA', external_neighbors=TRUE)
}
sapply(input_files, reducedDimNames)

umaps <- lapply(1:length(input_files), function(i){
  sce_umap <- reducedDim(x = input_files[[i]], type = "UMAP") %>%
  as.data.frame() %>%
  as_tibble() %>%
  bind_cols(colData(input_files[[i]]) %>% as_tibble()) %>%
  ggplot() +
  geom_point(aes(V1, V2, color=subsets_MT_percent)) +
  cowplot::theme_cowplot()+ggtitle(names(input_files)[i])
})
  
do.call('grid.arrange', umaps)

sce_denoise <- lapply(1:length(input_files), function(i){
  scran::denoisePCA(input_files[[i]], featuresec[[i]], subset.row=TopHVGs[[i]])
})
sce_denoise[[1]]

sce_denoise_umap <-  lapply(1:length(input_files), function(i){
  reducedDim(x = sce_denoise[[i]], type = "UMAP") %>%
  as.data.frame() %>%
  as_tibble() %>%
  bind_cols(colData(sce_denoise[[i]]) %>% as_tibble()) %>%
  ggplot() +
  geom_point(aes(V1, V2, color=subsets_MT_percent)) +
  cowplot::theme_cowplot()
})
do.call('grid.arrange', sce_denoise_umap)

umaps[[1]]
sce_denoise_umap[[1]]

#--------------------------------------------------------------------------------#

## doublets
dblts <- lapply(input_files, scDblFinder::scDblFinder)
dblts_plots <- lapply(1:length(input_files), function(i){
  ggplot(cbind.data.frame(reducedDim(x = input_files[[i]], type = "UMAP"),
                 col=dblts[[i]]$scDblFinder.score),
       aes(x=`1`, y=`2`, col=col))+
  geom_point()+
  jcolors::scale_color_jcolors_contin()
})
do.call('grid.arrange', dblts_plots)

emptydroplets_plots <- lapply(1:length(input_files), function(i){
  ggplot(cbind.data.frame(reducedDim(x = input_files[[i]], type = "UMAP"),
                          col=emptydrops[[i]]$LogProb),
         aes(x=`1`, y=`2`, col=col))+
    geom_point()+
    jcolors::scale_color_jcolors_contin()
})
do.call('grid.arrange', dblts_plots)

do.call('grid.arrange', emptydroplets_plots)
