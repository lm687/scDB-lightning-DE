# ## analyse DEG from Charly
# ## using some of his code
# ## but now with more stringent filtering; removing cells with fewer counts
# 
# rm(list = ls())
# setwd("/t1-data/project/cncb/shared/proj002/analyses/2020/10kCells/")
# source("~/projects/lighting/2_analyses/helper_functions.R")
# 
# library(Seurat)
# library(ggplot2)
# library(cowplot)
# library(plotly)
# library(reshape2)
# theme_set(theme_cowplot())
# 
# ##---------------------------------------------------------------------------------#
# ## Biomarkers
# 
# biomarkers_list <- readRDS("/t1-data/project/sims-lab/lmorrill/robjects/lightning/biomarkers_list.RDS")
# 
# Read10X("./FEB2019_G1_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G1_rep1.data
# Read10X("./FEB2019_G1_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G1_rep2.data
# Read10X("./FEB2019_G2_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G2_rep1.data
# Read10X("./FEB2019_G2_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G2_rep2.data
# Read10X("./FEB2019_G3_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G3_rep1.data
# Read10X("./FEB2019_G3_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G3_rep2.data
# Read10X("./FEB2019_G4_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G4_rep1.data
# Read10X("./FEB2019_G4_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G4_rep2.data
# 
# FEB2019_G1_rep1 <- CreateSeuratObject(FEB2019_G1_rep1.data, project = "FEB2019_G1_rep1", min.cells = 3, min.features = 200)
# FEB2019_G1_rep2 <- CreateSeuratObject(FEB2019_G1_rep2.data, project = "FEB2019_G1_rep2", min.cells = 3, min.features = 200)
# FEB2019_G2_rep1 <- CreateSeuratObject(FEB2019_G2_rep1.data, project = "FEB2019_G2_rep1", min.cells = 3, min.features = 200)
# FEB2019_G2_rep2 <- CreateSeuratObject(FEB2019_G2_rep2.data, project = "FEB2019_G2_rep2", min.cells = 3, min.features = 200)
# FEB2019_G3_rep1 <- CreateSeuratObject(FEB2019_G3_rep1.data, project = "FEB2019_G3_rep1", min.cells = 3, min.features = 200)
# FEB2019_G3_rep2 <- CreateSeuratObject(FEB2019_G3_rep2.data, project = "FEB2019_G3_rep2", min.cells = 3, min.features = 200)
# FEB2019_G4_rep1 <- CreateSeuratObject(FEB2019_G4_rep1.data, project = "FEB2019_G4_rep1", min.cells = 3, min.features = 200)
# FEB2019_G4_rep2 <- CreateSeuratObject(FEB2019_G4_rep2.data, project = "FEB2019_G4_rep2", min.cells = 3, min.features = 200)
# 
# rm(list=ls(pattern=".data"))
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
# # p1.1 <- FeatureScatter(FEB2019_G1_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # p1.2 <- FeatureScatter(FEB2019_G1_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # p1.3 <- FeatureScatter(FEB2019_G2_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # p1.4 <- FeatureScatter(FEB2019_G2_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # p1.5 <- FeatureScatter(FEB2019_G3_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # p1.6 <- FeatureScatter(FEB2019_G3_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # p1.7 <- FeatureScatter(FEB2019_G4_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # p1.8 <- FeatureScatter(FEB2019_G4_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # # CombinePlots(plots = list(p1.1, p1.2, p1.3, p1.4, p1.5, p1.6, p1.7, p1.8))
# # 
# # p2.1 <- FeatureScatter(FEB2019_G1_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # p2.2 <- FeatureScatter(FEB2019_G1_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # p2.3 <- FeatureScatter(FEB2019_G2_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # p2.4 <- FeatureScatter(FEB2019_G2_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # p2.5 <- FeatureScatter(FEB2019_G3_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # p2.6 <- FeatureScatter(FEB2019_G3_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # p2.7 <- FeatureScatter(FEB2019_G4_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # p2.8 <- FeatureScatter(FEB2019_G4_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # CombinePlots(plots = list(p2.1, p2.2, p2.3, p2.4, p2.5, p2.6, p2.7, p2.8))
# 
# FEB2019_G1_rep1 <- subset(FEB2019_G1_rep1, subset = nCount_RNA < 20000 & percent.mt < 15)
# FEB2019_G1_rep2 <- subset(FEB2019_G1_rep2, subset = nCount_RNA < 20000 & percent.mt < 15)
# FEB2019_G2_rep1 <- subset(FEB2019_G2_rep1, subset = nCount_RNA < 20000 & percent.mt < 15)
# FEB2019_G2_rep2 <- subset(FEB2019_G2_rep2, subset = nCount_RNA < 20000 & percent.mt < 15)
# FEB2019_G3_rep1 <- subset(FEB2019_G3_rep1, subset = nCount_RNA < 20000 & percent.mt < 15)
# FEB2019_G3_rep2 <- subset(FEB2019_G3_rep2, subset = nCount_RNA < 20000 & percent.mt < 15)
# FEB2019_G4_rep1 <- subset(FEB2019_G4_rep1, subset = nCount_RNA < 20000 & percent.mt < 15)
# FEB2019_G4_rep2 <- subset(FEB2019_G4_rep2, subset = nCount_RNA < 20000 & percent.mt < 15)
# 
# ## we are not subsetting cells with very few reads - possibly, those are subsetted by <cellranger count> already
# 
# samples <- apply(expand.grid(paste0('G', 1:4), paste0('_rep', 1:2)), 1, paste0, collapse='')
# samples <- gtools::mixedsort(samples)
# 
# par(mfrow=c(1,3), ask=T)
# sapply(paste0('FEB2019_', samples), function(name_object){
#   plot(density(log(get(name_object)$nCount_RNA))) ## bimodality seems to indicate a population of cells with very few reads
#   plot(density(log(get(name_object)$nFeature_RNA)))
#   plot((log(get(name_object)$nFeature_RNA)), log(get(name_object)$nCount_RNA), col=alpha('black', 0.02), main=name_object)
# })
# 
# FEB2019_G1_rep1 <- subset(FEB2019_G1_rep1, subset = nCount_RNA > exp(6.7))
# FEB2019_G1_rep2 <- subset(FEB2019_G1_rep2, subset = nCount_RNA > exp(6.8))
# FEB2019_G2_rep1 <- subset(FEB2019_G2_rep1, subset = nCount_RNA > exp(6.9))
# FEB2019_G2_rep2 <- subset(FEB2019_G2_rep2, subset = nCount_RNA > exp(6.8))
# FEB2019_G3_rep1 <- subset(FEB2019_G3_rep1, subset = nCount_RNA > exp(6.7)) ## strange one
# FEB2019_G3_rep2 <- subset(FEB2019_G3_rep2, subset = nCount_RNA > exp(6.7))
# FEB2019_G4_rep1 <- subset(FEB2019_G4_rep1, subset = nCount_RNA > exp(6.5)) ## very weird one
# FEB2019_G4_rep2 <- subset(FEB2019_G4_rep2, subset = nCount_RNA > exp(6.7)) ##??
# 
# par(mfrow=c(4,3), ask=T)
# sapply(paste0('FEB2019_', samples), function(name_object){
#   plot(density(log(get(name_object)$nCount_RNA))) ## bimodality seems to indicate a population of cells with very few reads
#   plot(density(log(get(name_object)$nFeature_RNA)))
#   plot((log(get(name_object)$nFeature_RNA)), log(get(name_object)$nCount_RNA), col=alpha('black', 0.02), main=name_object)
# })
# 
# FEB2019_G1_rep1 <- NormalizeData(FEB2019_G1_rep1)
# FEB2019_G1_rep2 <- NormalizeData(FEB2019_G1_rep2)
# FEB2019_G2_rep1 <- NormalizeData(FEB2019_G2_rep1)
# FEB2019_G2_rep2 <- NormalizeData(FEB2019_G2_rep2)
# FEB2019_G3_rep1 <- NormalizeData(FEB2019_G3_rep1)
# FEB2019_G3_rep2 <- NormalizeData(FEB2019_G3_rep2)
# FEB2019_G4_rep1 <- NormalizeData(FEB2019_G4_rep1)
# FEB2019_G4_rep2 <- NormalizeData(FEB2019_G4_rep2)
# 
# FEB2019_G1_rep1 <- FindVariableFeatures(FEB2019_G1_rep1)
# FEB2019_G1_rep2 <- FindVariableFeatures(FEB2019_G1_rep2)
# FEB2019_G2_rep1 <- FindVariableFeatures(FEB2019_G2_rep1)
# FEB2019_G2_rep2 <- FindVariableFeatures(FEB2019_G2_rep2)
# FEB2019_G3_rep1 <- FindVariableFeatures(FEB2019_G3_rep1)
# FEB2019_G3_rep2 <- FindVariableFeatures(FEB2019_G3_rep2)
# FEB2019_G4_rep1 <- FindVariableFeatures(FEB2019_G4_rep1)
# FEB2019_G4_rep2 <- FindVariableFeatures(FEB2019_G4_rep2)
# 
# FEB2019.anchors <- FindIntegrationAnchors(object.list = list(FEB2019_G1_rep1, FEB2019_G1_rep2, FEB2019_G2_rep1, FEB2019_G2_rep2, FEB2019_G3_rep1, FEB2019_G3_rep2, FEB2019_G4_rep1, FEB2019_G4_rep2), dims = 1:60)
# FEB2019.combined <- IntegrateData(anchorset = FEB2019.anchors, dims = 1:60)
# 
# DefaultAssay(FEB2019.combined) <- "integrated"
# FEB2019.combined <- ScaleData(FEB2019.combined)
# FEB2019.combined <- RunPCA(FEB2019.combined, npcs = 60)
# FEB2019.combined <- RunUMAP(FEB2019.combined, reduction = "pca", dims = 1:60)
# FEB2019.combined <- FindNeighbors(FEB2019.combined, reduction = "pca", dims = 1:60)
# FEB2019.combined <- FindClusters(FEB2019.combined, resolution = 4)
# FEB2019.combined <- RunTSNE(FEB2019.combined, reduction = "pca", dims = 1:60)
# 
# FEB2019.combined$stim <- ""
# FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G1_rep1"] <- "G1"
# FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G1_rep2"] <- "G1"
# FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G2_rep1"] <- "G2"
# FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G2_rep2"] <- "G2"
# FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G3_rep1"] <- "G3"
# FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G3_rep2"] <- "G3"
# FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G4_rep1"] <- "G4"
# FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G4_rep2"] <- "G4"
# 
# DefaultAssay(FEB2019.combined) <- "RNA"
# # MARKERS.FEB2019.combined <- FindAllMarkers(FEB2019.combined, logfc.threshold = 0.8, assay = "RNA")
# 
# ALL.average.expression <- AverageExpression(FEB2019.combined, verbose = F)
# 
# VlnPlot(FEB2019.combined, "Dsk", group.by = "stim", split.by = "orig.ident", pt.size = 0.1, slot = "counts", assay = "RNA", y.max = 10)
# 
# PCAPlot(FEB2019.combined, group.by= 'stim') ## this looks very strange
# 
# UMAPPlot(FEB2019.combined, group.by='stim') ## extremely similar to the UMAP of less stringent filtering
# TSNEPlot(FEB2019.combined, group.by='stim')
# 
# UMAPPlot(FEB2019.combined, group.by='nCount_RNA')
# FeaturePlot(FEB2019.combined, features='nCount_RNA')
# FeaturePlot(FEB2019.combined, features='nFeature_RNA')
# 
# # for(rep_it in samples){
# #   assign(x = paste0('FEB2019_', rep_it),
# #          value = add_dimred(get(paste0('FEB2019_', rep_it))))
# #   # assign(paste0('FEB2019_', rep_it), `[<-`(get(paste0('FEB2019_', rep_it))@meta.data[,'pca1'], 1, Reductions(get(paste0('FEB2019_', rep_it)), 'pca')@cell.embeddings[ match(rownames(get(paste0('FEB2019_', rep_it))@meta.data), 
# #   #                                                                                                                                                                           ownames(Reductions(get(paste0('FEB2019_', rep_it)), 'pca')@cell.embeddings)), 1]))
# #   assign(paste0('FEB2019_', rep_it), get(paste0('FEB2019_', rep_it))@meta.data, 1, cbind.data.frame(get(paste0('FEB2019_', rep_it))@meta.data, pca1=Reductions(get(paste0('FEB2019_', rep_it)), 'pca')@cell.embeddings[ match(rownames(get(paste0('FEB2019_', rep_it))@meta.data), 
# #                                                                                                                                                                             rownames(Reductions(get(paste0('FEB2019_', rep_it)), 'pca')@cell.embeddings)), 1]))
# #   
# # }
# 
# UMAPPlot(FEB2019_G1_rep1)
# TSNEPlot(FEB2019_G1_rep1)
# PCAPlot(FEB2019_G1_rep1) ## this looks VERY STRANGE
# 
# ## -------
# ## I couldn't automate this (see above)
# FEB2019_G1_rep1@meta.data[,'pca1'] <- Reductions(FEB2019_G1_rep1, 'pca')@cell.embeddings[ match(rownames(FEB2019_G1_rep1@meta.data), rownames(Reductions(FEB2019_G1_rep1, 'pca')@cell.embeddings)), 1]
# FEB2019_G1_rep1@meta.data[,'pca2'] <- Reductions(FEB2019_G1_rep1, 'pca')@cell.embeddings[ match(rownames(FEB2019_G1_rep1@meta.data), rownames(Reductions(FEB2019_G1_rep1, 'pca')@cell.embeddings)), 2]
# FEB2019_G1_rep2@meta.data[,'pca1'] <- Reductions(FEB2019_G1_rep2, 'pca')@cell.embeddings[ match(rownames(FEB2019_G1_rep2@meta.data), rownames(Reductions(FEB2019_G1_rep2, 'pca')@cell.embeddings)), 1]
# FEB2019_G1_rep2@meta.data[,'pca2'] <- Reductions(FEB2019_G1_rep2, 'pca')@cell.embeddings[ match(rownames(FEB2019_G1_rep2@meta.data), rownames(Reductions(FEB2019_G1_rep2, 'pca')@cell.embeddings)), 2]
# 
# FEB2019_G2_rep1@meta.data[,'pca1'] <- Reductions(FEB2019_G2_rep1, 'pca')@cell.embeddings[ match(rownames(FEB2019_G2_rep1@meta.data), rownames(Reductions(FEB2019_G2_rep1, 'pca')@cell.embeddings)), 1]
# FEB2019_G2_rep1@meta.data[,'pca2'] <- Reductions(FEB2019_G2_rep1, 'pca')@cell.embeddings[ match(rownames(FEB2019_G2_rep1@meta.data), rownames(Reductions(FEB2019_G2_rep1, 'pca')@cell.embeddings)), 2]
# FEB2019_G2_rep2@meta.data[,'pca1'] <- Reductions(FEB2019_G2_rep2, 'pca')@cell.embeddings[ match(rownames(FEB2019_G2_rep2@meta.data), rownames(Reductions(FEB2019_G2_rep2, 'pca')@cell.embeddings)), 1]
# FEB2019_G2_rep2@meta.data[,'pca2'] <- Reductions(FEB2019_G2_rep2, 'pca')@cell.embeddings[ match(rownames(FEB2019_G2_rep2@meta.data), rownames(Reductions(FEB2019_G2_rep2, 'pca')@cell.embeddings)), 2]
# 
# FEB2019_G3_rep1@meta.data[,'pca1'] <- Reductions(FEB2019_G3_rep1, 'pca')@cell.embeddings[ match(rownames(FEB2019_G3_rep1@meta.data), rownames(Reductions(FEB2019_G3_rep1, 'pca')@cell.embeddings)), 1]
# FEB2019_G3_rep1@meta.data[,'pca2'] <- Reductions(FEB2019_G3_rep1, 'pca')@cell.embeddings[ match(rownames(FEB2019_G3_rep1@meta.data), rownames(Reductions(FEB2019_G3_rep1, 'pca')@cell.embeddings)), 2]
# FEB2019_G3_rep2@meta.data[,'pca1'] <- Reductions(FEB2019_G3_rep2, 'pca')@cell.embeddings[ match(rownames(FEB2019_G3_rep2@meta.data), rownames(Reductions(FEB2019_G3_rep2, 'pca')@cell.embeddings)), 1]
# FEB2019_G3_rep2@meta.data[,'pca2'] <- Reductions(FEB2019_G3_rep2, 'pca')@cell.embeddings[ match(rownames(FEB2019_G3_rep2@meta.data), rownames(Reductions(FEB2019_G3_rep2, 'pca')@cell.embeddings)), 2]
# 
# FEB2019_G4_rep1@meta.data[,'pca1'] <- Reductions(FEB2019_G4_rep1, 'pca')@cell.embeddings[ match(rownames(FEB2019_G4_rep1@meta.data), rownames(Reductions(FEB2019_G4_rep1, 'pca')@cell.embeddings)), 1]
# FEB2019_G4_rep1@meta.data[,'pca2'] <- Reductions(FEB2019_G4_rep1, 'pca')@cell.embeddings[ match(rownames(FEB2019_G4_rep1@meta.data), rownames(Reductions(FEB2019_G4_rep1, 'pca')@cell.embeddings)), 2]
# FEB2019_G4_rep2@meta.data[,'pca1'] <- Reductions(FEB2019_G4_rep2, 'pca')@cell.embeddings[ match(rownames(FEB2019_G4_rep2@meta.data), rownames(Reductions(FEB2019_G4_rep2, 'pca')@cell.embeddings)), 1]
# FEB2019_G4_rep2@meta.data[,'pca2'] <- Reductions(FEB2019_G4_rep2, 'pca')@cell.embeddings[ match(rownames(FEB2019_G4_rep2@meta.data), rownames(Reductions(FEB2019_G4_rep2, 'pca')@cell.embeddings)), 2]
# 
# ## -------
# 
# library(gridExtra)
# do.call('grid.arrange', lapply(samples, function(i) FeaturePlot(get(paste0('FEB2019_', i)), features='pca1')+labs(title = i)))
# do.call('grid.arrange', lapply(samples, function(i) FeaturePlot(get(paste0('FEB2019_', i)), features='pca1')+labs(title = i)))
# do.call('grid.arrange', lapply(samples, function(i) FeaturePlot(get(paste0('FEB2019_', i)), reduction = 'pca', features = 'nCount_RNA')+labs(title = i)))
# par(mfrow=c(4,2))
# lapply(samples, function(i) plot(get(paste0('FEB2019_', i))$nCount_RNA,
#                                                          get(paste0('FEB2019_', i))@meta.data[,'pca1']))
# ## pca is not done on normalised counts? or the normalisation is not appropriate
# 
# FeaturePlot(FEB2019_G1_rep1, features = Reductions(FEB2019_G1_rep1, 'pca')[,1])
# 
# FeaturePlot(FEB2019_G1_rep1, features='lncRNA:roX1')
# 
# ## -------
# ## normalise using some other method
# 
# for(rep_it in samples){
#   assign(x = paste0('FEB2019_', rep_it), SCTransform(get(paste0('FEB2019_', rep_it))))
# }
# 
# ## re-compute dim red with SCT normalisation
# for(rep_it in samples){
#   assign(x = paste0('FEB2019_', rep_it), add_dimred(get(paste0('FEB2019_', rep_it))))
# }
# 
# do.call('grid.arrange', lapply(samples, function(i) FeaturePlot(get(paste0('FEB2019_', i)), features='nCount_RNA')+labs(title = i)))
# do.call('grid.arrange', lapply(samples, function(i) FeaturePlot(get(paste0('FEB2019_', i)), features='nCount_RNA', reduction = 'pca')+labs(title = i)))
# do.call('grid.arrange', lapply(samples, function(i) FeaturePlot(get(paste0('FEB2019_', i)), features='nCount_RNA', reduction = 'tsne')+labs(title = i)))
# 
# FeaturePlot(FEB2019_G1_rep1, features = Reductions(FEB2019_G1_rep1, 'pca')[,1])
# PCAPlot(FEB2019_G1_rep1)
# UMAPPlot(FEB2019_G1_rep1) ## horrible
# TSNEPlot(FEB2019_G1_rep1)
# 
# FeaturePlot(FEB2019_G1_rep1, features='nCount_RNA')
# FeaturePlot(FEB2019_G1_rep1, features='nCount_RNA', reduction = 'pca')
# 
# FEB2019_G1_rep1 <- RunPCA(FEB2019_G1_rep1, npcs = 60, assay = "SCT" )
# FEB2019_G1_rep1
# PCAPlot(FEB2019_G1_rep1)
# FeaturePlot(FEB2019_G1_rep1, features = "nCount_RNA") ## no relationship between pc and rna count
# 
