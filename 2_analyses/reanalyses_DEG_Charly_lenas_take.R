## analyse DEG from Charly (using his cellranger counts)
## and using some of his code
## with removal of empty drops and doublets, and other changes made by Lena

rm(list = ls())
setwd("/t1-data/project/cncb/shared/proj002/analyses/2020/10kCells/")

library(Seurat)
library(ggplot2)
library(cowplot)
library(plotly)
library(reshape2)
library(gridExtra)
library(DropletUtils)
theme_set(theme_cowplot())

samples <- apply(expand.grid(paste0('G', 1:4), paste0('_rep', 1:2)), 1, paste0, collapse='')
samples

for(rep_it in samples){
  assign(x = paste0('FEB2019_', rep_it, '.data'),
         value = Read10X(paste0('./FEB2019_', rep_it, '/outs/filtered_feature_bc_matrix')))
}

for(rep_it in samples){
  assign(x = paste0('FEB2019_', rep_it),
         value = CreateSeuratObject(get(paste0('FEB2019_', rep_it, '.data')), project = paste0('FEB2019_', rep_it),
                                    min.cells = 3, min.features = 200))
}

rm(list=ls(pattern=".data"))

# COPYFEB2019_G1_rep1 <- FEB2019_G1_rep1
# rep_it= samples[1]
# for(rep_it in samples){
#   assign(x = paste0('COPYFEB2019_', rep_it), pos = 'percent.mt',
#          PercentageFeatureSet(get(paste0('FEB2019_', rep_it)), pattern = "^mt:"))
#   # assign(x = paste0('COPYFEB2019_', rep_it, "[['percent.mt']]"), 
#   #        PercentageFeatureSet(get(paste0('FEB2019_', rep_it)), pattern = "^mt:"))
#   assign(x = paste0('COPYFEB2019_', rep_it, "@meta.data"), 
#          cbind.data.frame(get(paste0('COPYFEB2019_', rep_it))@meta.data,
#                           percent.mt=PercentageFeatureSet(get(paste0('FEB2019_', rep_it)), pattern = "^mt:")[,1])) ## not enough
#   get(paste0('COPYFEB2019_', rep_it))[['percent.mt']] <- PercentageFeatureSet(get(paste0('FEB2019_', rep_it)), pattern = "^mt:") ## doesn't work
#   # assign(x = paste0('FEB2019_', rep_it), 
#   #        list(paste0('FEB2019_', rep_it, collapse=''), percent.mt=PercentageFeatureSet(get(paste0('FEB2019_', rep_it)), pattern = "^mt:")))
# }

FEB2019_G1_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G1_rep1, pattern = "^mt:")
FEB2019_G1_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G1_rep2, pattern = "^mt:")
FEB2019_G2_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G2_rep1, pattern = "^mt:")
FEB2019_G2_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G2_rep2, pattern = "^mt:")
FEB2019_G3_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G3_rep1, pattern = "^mt:")
FEB2019_G3_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G3_rep2, pattern = "^mt:")
FEB2019_G4_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G4_rep1, pattern = "^mt:")
FEB2019_G4_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G4_rep2, pattern = "^mt:")

do.call('grid.arrange', lapply(samples, function(rep_it) FeatureScatter(get(paste0('FEB2019_', rep_it)),
                                                                        feature1 = "nCount_RNA", feature2 = "percent.mt")))

do.call('grid.arrange', lapply(samples, function(rep_it) FeatureScatter(get(paste0('FEB2019_', rep_it)),
                                                                        feature1 = "nCount_RNA", feature2 = "nFeature_RNA")))

for(rep_it in samples){
  assign(x = paste0('FEB2019_', rep_it, '_summarisedexp'),
         value = DropletUtils::read10xCounts(paste0('./FEB2019_', rep_it, '/outs/raw_feature_bc_matrix')))
}
## Empty droplets, doublets
emptydrops <- lapply(samples, function(rep_it) DropletUtils::emptyDrops(m = as(assay(FEB2019_G1_rep1_summarisedexp), 'matrix')))
DropletUtils::emptyDrops(m = FEB2019_G1_rep1_summarisedexp)

## Filter by other means
for(rep_it in samples){
  assign(x = paste0('FEB2019_', rep_it), subset(get(paste0('FEB2019_', rep_it)), subset = nCount_RNA < 20000 & percent.mt < 15))
}

