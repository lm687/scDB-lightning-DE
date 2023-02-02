## analyse DEG from Charly
## using some of his code

rm(list = ls())
setwd("/t1-data/project/cncb/shared/proj002/analyses/2020/10kCells/")

library(Seurat)
library(ggplot2)
library(cowplot)
library(plotly)
library(reshape2)
theme_set(theme_cowplot())

source("~/projects/lighting/2_analyses/helper_functions.R")

##---------------------------------------------------------------------------------#
## Biomarkers

biomarkers_list <- readRDS("/t1-data/project/sims-lab/lmorrill/robjects/lightning/biomarkers_list.RDS")

Read10X("./FEB2019_G1_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G1_rep1.data
Read10X("./FEB2019_G1_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G1_rep2.data
Read10X("./FEB2019_G2_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G2_rep1.data
Read10X("./FEB2019_G2_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G2_rep2.data
Read10X("./FEB2019_G3_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G3_rep1.data
Read10X("./FEB2019_G3_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G3_rep2.data
Read10X("./FEB2019_G4_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G4_rep1.data
Read10X("./FEB2019_G4_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G4_rep2.data

FEB2019_G1_rep1 <- CreateSeuratObject(FEB2019_G1_rep1.data, project = "FEB2019_G1_rep1", min.cells = 3, min.features = 200)
FEB2019_G1_rep2 <- CreateSeuratObject(FEB2019_G1_rep2.data, project = "FEB2019_G1_rep2", min.cells = 3, min.features = 200)
FEB2019_G2_rep1 <- CreateSeuratObject(FEB2019_G2_rep1.data, project = "FEB2019_G2_rep1", min.cells = 3, min.features = 200)
FEB2019_G2_rep2 <- CreateSeuratObject(FEB2019_G2_rep2.data, project = "FEB2019_G2_rep2", min.cells = 3, min.features = 200)
FEB2019_G3_rep1 <- CreateSeuratObject(FEB2019_G3_rep1.data, project = "FEB2019_G3_rep1", min.cells = 3, min.features = 200)
FEB2019_G3_rep2 <- CreateSeuratObject(FEB2019_G3_rep2.data, project = "FEB2019_G3_rep2", min.cells = 3, min.features = 200)
FEB2019_G4_rep1 <- CreateSeuratObject(FEB2019_G4_rep1.data, project = "FEB2019_G4_rep1", min.cells = 3, min.features = 200)
FEB2019_G4_rep2 <- CreateSeuratObject(FEB2019_G4_rep2.data, project = "FEB2019_G4_rep2", min.cells = 3, min.features = 200)

rm(list=ls(pattern=".data"))

FEB2019_G1_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G1_rep1, pattern = "^mt:")
FEB2019_G1_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G1_rep2, pattern = "^mt:")
FEB2019_G2_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G2_rep1, pattern = "^mt:")
FEB2019_G2_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G2_rep2, pattern = "^mt:")
FEB2019_G3_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G3_rep1, pattern = "^mt:")
FEB2019_G3_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G3_rep2, pattern = "^mt:")
FEB2019_G4_rep1[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G4_rep1, pattern = "^mt:")
FEB2019_G4_rep2[["percent.mt"]] <- PercentageFeatureSet(FEB2019_G4_rep2, pattern = "^mt:")

p1.1 <- FeatureScatter(FEB2019_G1_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
p1.2 <- FeatureScatter(FEB2019_G1_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p1.3 <- FeatureScatter(FEB2019_G2_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
p1.4 <- FeatureScatter(FEB2019_G2_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p1.5 <- FeatureScatter(FEB2019_G3_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
p1.6 <- FeatureScatter(FEB2019_G3_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt")
p1.7 <- FeatureScatter(FEB2019_G4_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
p1.8 <- FeatureScatter(FEB2019_G4_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt")
CombinePlots(plots = list(p1.1, p1.2, p1.3, p1.4, p1.5, p1.6, p1.7, p1.8))

p2.1 <- FeatureScatter(FEB2019_G1_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2.2 <- FeatureScatter(FEB2019_G1_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2.3 <- FeatureScatter(FEB2019_G2_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2.4 <- FeatureScatter(FEB2019_G2_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2.5 <- FeatureScatter(FEB2019_G3_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2.6 <- FeatureScatter(FEB2019_G3_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2.7 <- FeatureScatter(FEB2019_G4_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2.8 <- FeatureScatter(FEB2019_G4_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(p2.1, p2.2, p2.3, p2.4, p2.5, p2.6, p2.7, p2.8))

FEB2019_G1_rep1 <- subset(FEB2019_G1_rep1, subset = nCount_RNA < 20000 & percent.mt < 15)
FEB2019_G1_rep2 <- subset(FEB2019_G1_rep2, subset = nCount_RNA < 20000 & percent.mt < 15)
FEB2019_G2_rep1 <- subset(FEB2019_G2_rep1, subset = nCount_RNA < 20000 & percent.mt < 15)
FEB2019_G2_rep2 <- subset(FEB2019_G2_rep2, subset = nCount_RNA < 20000 & percent.mt < 15)
FEB2019_G3_rep1 <- subset(FEB2019_G3_rep1, subset = nCount_RNA < 20000 & percent.mt < 15)
FEB2019_G3_rep2 <- subset(FEB2019_G3_rep2, subset = nCount_RNA < 20000 & percent.mt < 15)
FEB2019_G4_rep1 <- subset(FEB2019_G4_rep1, subset = nCount_RNA < 20000 & percent.mt < 15)
FEB2019_G4_rep2 <- subset(FEB2019_G4_rep2, subset = nCount_RNA < 20000 & percent.mt < 15)

## we are not subsetting cells with very few reads - possibly, those are subsetted by <cellranger count> already

plot(density(FEB2019_G1_rep1$nCount_RNA))
plot(density(log(FEB2019_G1_rep1$nCount_RNA))) ## bimodality seems to indicate a population of cells with very few reads
plot(density(log(FEB2019_G1_rep1$nFeature_RNA)))
plot((log(FEB2019_G1_rep1$nFeature_RNA)), log(FEB2019_G1_rep1$nCount_RNA), col=alpha('black', 0.02))


FEB2019_G1_rep1 <- NormalizeData(FEB2019_G1_rep1)
FEB2019_G1_rep2 <- NormalizeData(FEB2019_G1_rep2)
FEB2019_G2_rep1 <- NormalizeData(FEB2019_G2_rep1)
FEB2019_G2_rep2 <- NormalizeData(FEB2019_G2_rep2)
FEB2019_G3_rep1 <- NormalizeData(FEB2019_G3_rep1)
FEB2019_G3_rep2 <- NormalizeData(FEB2019_G3_rep2)
FEB2019_G4_rep1 <- NormalizeData(FEB2019_G4_rep1)
FEB2019_G4_rep2 <- NormalizeData(FEB2019_G4_rep2)

FEB2019_G1_rep1 <- FindVariableFeatures(FEB2019_G1_rep1)
FEB2019_G1_rep2 <- FindVariableFeatures(FEB2019_G1_rep2)
FEB2019_G2_rep1 <- FindVariableFeatures(FEB2019_G2_rep1)
FEB2019_G2_rep2 <- FindVariableFeatures(FEB2019_G2_rep2)
FEB2019_G3_rep1 <- FindVariableFeatures(FEB2019_G3_rep1)
FEB2019_G3_rep2 <- FindVariableFeatures(FEB2019_G3_rep2)
FEB2019_G4_rep1 <- FindVariableFeatures(FEB2019_G4_rep1)
FEB2019_G4_rep2 <- FindVariableFeatures(FEB2019_G4_rep2)

FEB2019.anchors <- FindIntegrationAnchors(object.list = list(FEB2019_G1_rep1, FEB2019_G1_rep2, FEB2019_G2_rep1, FEB2019_G2_rep2, FEB2019_G3_rep1, FEB2019_G3_rep2, FEB2019_G4_rep1, FEB2019_G4_rep2), dims = 1:60)
FEB2019.combined <- IntegrateData(anchorset = FEB2019.anchors, dims = 1:60)

DefaultAssay(FEB2019.combined) <- "integrated"
FEB2019.combined <- ScaleData(FEB2019.combined)
FEB2019.combined <- RunPCA(FEB2019.combined, npcs = 60)
FEB2019.combined <- RunUMAP(FEB2019.combined, reduction = "pca", dims = 1:60)
FEB2019.combined <- FindNeighbors(FEB2019.combined, reduction = "pca", dims = 1:60)
FEB2019.combined <- FindClusters(FEB2019.combined, resolution = 4)
FEB2019.combined <- RunTSNE(FEB2019.combined, reduction = "pca", dims = 1:60)

FEB2019.combined$stim <- ""
FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G1_rep1"] <- "G1"
FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G1_rep2"] <- "G1"
FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G2_rep1"] <- "G2"
FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G2_rep2"] <- "G2"
FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G3_rep1"] <- "G3"
FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G3_rep2"] <- "G3"
FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G4_rep1"] <- "G4"
FEB2019.combined$stim[FEB2019.combined[["orig.ident"]] == "FEB2019_G4_rep2"] <- "G4"

DefaultAssay(FEB2019.combined) <- "RNA"
MARKERS.FEB2019.combined <- FindAllMarkers(FEB2019.combined, logfc.threshold = 0.8, assay = "RNA")
save(MARKERS.FEB2019.combined, file = "./FEB2019.MARKERS.allgroups.Robj")


ALL.average.expression <- AverageExpression(FEB2019.combined, verbose = F)

grep("^TE-*", rownames(FEB2019.combined@assays$RNA), value = T) -> telist ## LM: PROBLEM. telist is created after it was called below. hence this line has been moved
ALL.avg.TEexpression <- log1p(ALL.average.expression$RNA[telist,]) ### PROBLEM. telist is needed
ALL.avg.TEexpression$gene <- rownames(ALL.avg.TEexpression)
p1 <- ggplot(ALL.avg.TEexpression, aes(G1, G2)) + geom_point() + ggtitle("G1 vs. G2") ## ERROR

Idents(FEB2019.combined) <- FEB2019.combined$integrated_snn_res.4
KCs.cell <- subset(FEB2019.combined, idents = c("3", "7", "29", "32", "36", "78", "90", "100", "112"))
Idents(KCs.cell) <- "stim"
avg.KCs.cell <- log1p(AverageExpression(KCs.cell, verbose = F)$RNA)
avg.KCs.cell$gene <- rownames(avg.KCs.cell)
p1 <- ggplot(avg.KCs.cell, aes(G3, G2, name = gene)) + geom_point() + ggtitle("G3 vs. G2")
ggplotly(p1)

p2 <- ggplot(avg.KCs.cell[telist,], aes(G3, G2, name = gene)) + geom_point() + ggtitle("G3 vs. G2")


VlnPlot(FEB2019.combined, "Dsk", group.by = "stim", split.by = "orig.ident", pt.size = 0.1, slot = "counts", assay = "RNA", y.max = 10)


for (i in telist) {
  jpeg(file = paste("~/Dropbox/CloudDesktop/TEplots/TE_count_", i, ".jpeg", sep=""))
  plot1 <- VlnPlot(FEB2019.combined, idents = c(27), group.by = "orig.ident",features = as.character(i), slot = "counts", assay = "RNA")
  plot(plot1)
  dev.off()
}


save.image(file = "~/projects/lighting/data/robjects/image_DEA_Charly.RData")

## --------------------------------------------------------- ##
## FROM NOW ON: Lena
## --------------------------------------------------------- ##
# load(file = "~/projects/lighting/data/robjects/image_DEA_Charly.RData")

## --------------------------------------------------------- ##
## what are the Chrimson positive cells?
#' es, absolutely: It’s mainly a bunch of positively reinforcing dopamine neurons, plus a group of off-target 
#' cholinergic neurons (plus cells in the eye, but they are not part of the data set) . The driver we used is 0273-Gal4

FEB2019.combined['Chrimson']
plotdensity <- function(arg_vec, log2p1lusimp=T, value_imput=1){
  x <- (reshape2::melt(arg_vec))
  if(log2p1lusimp){
    ggplot(x, aes(x=value+value_imput))+geom_density()+scale_x_continuous(trans = "log10")
  }else{
    ggplot(x, aes(x=value))+geom_density()+scale_x_continuous(trans = "log2")
  }
}
plot(density(as.vector(FEB2019.combined@assays$RNA['Chrimson'])))
plotdensity(as.vector(FEB2019.combined@assays$RNA['Chrimson']), value_imput = 0.0000001)
plotdensity(as.vector(FEB2019.combined@assays$RNA['Chrimson']), value_imput = 0)
FEB2019.combined@assays$integrated['Chrimson']
max(as.vector(FEB2019.combined@assays$RNA['Chrimson']))
plotdensity(as.vector(FEB2019.combined@assays$RNA['Chrimson']), value_imput = 0)


## --------------------------------------------------------- ##
## umaps

UMAPPlot(FEB2019.combined, group.by='stim')
ggsave("~/projects/lighting/3_results/stim_umap.png", height = 4, width = 4)

splitUMAPPlot(FEB2019.combined, group.by='stim')
ggsave("~/projects/lighting/3_results/stim_umap_split.png", height = 6, width = 8)

FeaturePlot(FEB2019.combined, features = "Chrimson")
ggsave("~/projects/lighting/3_results/chrimsom_umap.png", height = 4, width = 4)

FeaturePlot(FEB2019.combined, features = "Chrimson", reduction = "tsne")
ggsave("~/projects/lighting/3_results/chrimsom_tsne.png", height = 4, width = 4)

DimPlot(FEB2019.combined, reduction = "umap")+ NoLegend()
ggsave("~/projects/lighting/3_results/clusters_umap.png", height = 4, width = 4)

DimPlot(FEB2019.combined, reduction = "tsne")+ NoLegend()
ggsave("~/projects/lighting/3_results/clusters_tsne.png", height = 4, width = 4)

FEB2019.combined@assays

FeaturePlot(FEB2019.combined, features = biomarkers_list$`gene name`)

VlnPlot(FEB2019.combined, features  = biomarkers_list$`gene name`)

# aggregate(x = as(FEB2019.combined@assays$RNA@counts['Vmat',], 'vector'), by = lapply(unique( FEB2019.combined$seurat_clusters), function(i)  FEB2019.combined$seurat_clusters == i ), FUN='mean')
Vmat_across_clusters <- aggregate(x = as(FEB2019.combined@assays$RNA@counts['Vmat',], 'vector'), by = list(clusters=FEB2019.combined$seurat_clusters), FUN='mean')
ggplot(Vmat_across_clusters, aes(x=factor(clusters, levels=clusters[order(x)]), y=x))+geom_point()+geom_line()

biomarkers_list_vec <- c(biomarkers_list$`gene name`, 'Vmat', 'DAT', 'Frq1', 'otp', 'Lmx1a', 'Tdc2', 'Tbh', 'Tbh', 'kek1', 'hth', 'mirr', 'vvl', 'ey', 'Lim3')
## aggregate by mean
# biomarkers_across_clusters <- sapply(biomarkers_list_vec, function(biomarker_it){
#   aggregate(x = as(FEB2019.combined@assays$RNA@counts[biomarker_it,], 'vector'), by = list(clusters=FEB2019.combined$seurat_clusters), FUN='mean')[,2]
# })
## aggregate by median
biomarkers_across_clusters <- sapply(biomarkers_list_vec, function(biomarker_it){
  aggregate(x = as(FEB2019.combined@assays$RNA@counts[biomarker_it,], 'vector'), by = list(clusters=FEB2019.combined$seurat_clusters), FUN='median')[,2]
})
pheatmap::pheatmap(biomarkers_across_clusters)

# Seurat::clust (FEB2019.combined)

names_samples <- c('FEB2019_G1_rep1',
                   'FEB2019_G1_rep2',
                   'FEB2019_G2_rep1',
                   'FEB2019_G2_rep2',
                   'FEB2019_G3_rep1',
                   'FEB2019_G3_rep2',
                   'FEB2019_G4_rep1',
                   'FEB2019_G4_rep2')
dimred_individual <- lapply(c(FEB2019_G1_rep1,
                                   FEB2019_G1_rep2,
                                   FEB2019_G2_rep1,
                                   FEB2019_G2_rep2,
                                   FEB2019_G3_rep1,
                                   FEB2019_G3_rep2,
                                   FEB2019_G4_rep1,
                                   FEB2019_G4_rep2), add_dimred)
names(dimred_individual) <- names_samples
for(i in names(dimred_individual)){
  FeaturePlot(dimred_individual[[i]], features = "Chrimson")
  ggsave(paste0("~/projects/lighting/3_results/chrimsom_umap_", i, ".png"), height = 4, width = 4)
}

rownames(FEB2019.combined$nCount_RNA)

Seurat::Assays(FEB2019.combined)

FEB2019.combined@assays$RNA@counts['Chrimson',]
plot(density(FEB2019.combined@assays$RNA@counts['Chrimson',]))
plot(density(log(FEB2019.combined@assays$RNA@counts['Chrimson',])))
plot(density(log(1+FEB2019.combined@assays$RNA@counts['Chrimson',])))
plot(density(as.vector(FetchData(FEB2019.combined, 'Chrimson'))[[1]]))
table(as.vector(FetchData(FEB2019.combined, 'Chrimson'))[[1]] > 0)

## hard threshold of at least one count in Chrimson (i.e. probably many false negatives)
table(FEB2019.combined@assays$RNA@counts['Chrimson',] > 0)

##---------------------------------------------------------------------------------#
## Split by Chrimson

Chrimsonpos <- FEB2019.combined@assays$RNA@counts['Chrimson',] > 0
table(Chrimsonpos)

Chrimsonposcells <- subset(FEB2019.combined, cells = which(Chrimsonpos))
Chrimsonnegcells <- subset(FEB2019.combined, cells = which(!Chrimsonpos))
Chrimsonposcells ## 1899 samples
Chrimsonposcells <- FindVariableFeatures(Chrimsonposcells)
Chrimsonposcells <- add_dimred(Chrimsonposcells)
Chrimsonnegcells <- FindVariableFeatures(Chrimsonnegcells)
Chrimsonnegcells <- add_dimred(Chrimsonnegcells)
saveRDS(Chrimsonposcells, "~/projects/lighting/data/robjects/Chrimsonposcells.RDS")
saveRDS(Chrimsonnegcells, "~/projects/lighting/data/robjects/Chrimsonnegcells.RDS")

# Chrimsonposcells <- readRDS("~/projects/lighting/data/robjects/Chrimsonposcells.RDS")
# Chrimsonnegcells <- readRDS("~/projects/lighting/data/robjects/Chrimsonnegcells.RDS")

##---------------------------------------------------------------------------------#
## Biomarkers
## top markers for each of the clusters
xtable::xtable(MARKERS.FEB2019.combined[sapply(unique(MARKERS.FEB2019.combined$cluster), function(i) which(MARKERS.FEB2019.combined$cluster == i)[1]),])

UMAPPlot(Chrimsonposcells, group.by='stim')
UMAPPlot(Chrimsonposcells)
FeaturePlot(Chrimsonposcells, features='Chrimson')
UMAPPlot(Chrimsonnegcells, group.by='stim')
UMAPPlot(Chrimsonnegcells)
TSNEPlot(Chrimsonposcells)

FeaturePlot(Chrimsonposcells, features='nCount_RNA')
FeaturePlot(Chrimsonposcells, features='nFeature_RNA')

plot(Chrimsonposcells$nFeature_RNA, Chrimsonposcells$nCount_RNA)

my_VlnPlot <- function(seurat_object, features, aggregation_fun='median', arg_ncol=4, logtrans=T){
  # biomarkers_across_clusters <- sapply(features, function(biomarker_it){
  #   aggregate(x = as(seurat_object@assays$RNA@counts[biomarker_it,], 'vector'), by = list(clusters=seurat_object$seurat_clusters), FUN=aggregation_fun)[,2]
  # })
  # rownames(biomarkers_across_clusters) <- levels(seurat_object$seurat_clusters)
  subsetbiomarkers <- reshape2::melt(as(seurat_object@assays$RNA@counts[features,], 'matrix'))
  subsetbiomarkers$cluster <- seurat_object$seurat_clusters[match(subsetbiomarkers$Var2, colnames(seurat_object@assays$RNA@counts))]
  a <- ggplot(subsetbiomarkers, aes(x=cluster, y=value, group=cluster))+
    geom_boxplot()+
    geom_jitter(alpha=0.2)+
    facet_wrap(.~Var1, ncol=arg_ncol, scales = "free_y")
  if(logtrans){
    a <- a + scale_y_continuous(trans = "log2")
  }
  a
}

Chrimsonposcells_highlevel_cluster <- FindClusters(Chrimsonposcells, resolution = .1)
UMAPPlot(Chrimsonposcells_highlevel_cluster)
VlnPlot(Chrimsonposcells_highlevel_cluster, features = biomarkers_list_vec[1:5])
my_VlnPlot(Chrimsonposcells_highlevel_cluster, features = biomarkers_list_vec[1:5])
my_VlnPlot(Chrimsonposcells_highlevel_cluster, features = biomarkers_list_vec)
my_VlnPlot(Chrimsonposcells_highlevel_cluster, features = biomarkers_list_vec, logtrans = F, arg_ncol = 7)

for(i in 1:nrow(biomarkers_list)){
  cat(biomarkers_list[i,'gene name'], '\n')
  FeaturePlot(Chrimsonposcells, features=biomarkers_list[i,'gene name'])+labs(title = paste0(biomarkers_list[i,'gene name'],
                                                                                             "\n",
                                                                                             biomarkers_list[i,'marker']))
  ggsave(paste0("~/projects/lighting/3_results/markers/",biomarkers_list[i,'gene name'], ".png"), height = 4, width = 4)
  FeaturePlot(Chrimsonposcells, features=biomarkers_list[i,'gene name'], reduction = "tsne")+
    labs(title = paste0(biomarkers_list[i,'gene name'],
                                                                                             "\n",
                                                                                             biomarkers_list[i,'marker']))
  ggsave(paste0("~/projects/lighting/3_results/markers/TSNE_",biomarkers_list[i,'gene name'], ".png"), height = 4, width = 4)
  
}

normalise_rw <- function(i){
  sweep(i, 1, rowSums(i), '/')
}
rownames_to_col <- function(i){
  data.frame(names=rownames(i), i)
}

relevel_by_value_column <- function(i){
  i$names <- factor(i$names, levels=i$names[order(i$value)])
  return(i)
}

## which cell clusters are in each of the two (artificial groups) according to the UMAP?
# ggplot(relevel_by_value_column(rownames_to_col(reshape2::melt(sort(normalise_rw(table(Chrimsonposcells@meta.data$seurat_clusters,
#       Chrimsonposcells@meta.data$seurat_clusters_two_groups))[,1])))),
#       aes(x=names, y=1-value, group=1))+geom_point()+geom_line()+theme_bw()+
#   labs(x='Cluster name', y='Fraction of Chrimson+ cells')
# ggsave("~/projects/lighting/3_results/chrimsonpos/fractioncells_clusters_two_groups_cluster.png", height = 2, width = 7)

## genes characterising each cluster
Seurat::DimHeatmap(Chrimsonposcells)
Seurat::DoHeatmap(Chrimsonposcells)
## DE between two cell groups
Seurat::FoldChange(Chrimsonposcells, ident.1='1', ident.2='2')
give_top_logf_genes(Seurat::FoldChange(Chrimsonposcells, ident.1='1', ident.2='3'))
Seurat::FoldChange(Chrimsonposcells, ident.1='1', ident.2='4')
AllMarkersChrimsonPos <- Seurat::FindAllMarkers(Chrimsonposcells)
AllMarkersChrimsonPos
xtable::xtable(AllMarkersChrimsonPos[sapply(unique(AllMarkersChrimsonPos$cluster), function(i) which(AllMarkersChrimsonPos$cluster == i)[1]),])

Seurat::TopFeatures(Chrimsonposcells)

Chrimsonposcells@assays$RNA@meta.features
Chrimsonposcells@meta.data$seurat_clusters_two_groups <- Chrimsonposcells@reductions$umap@cell.embeddings[,1] < -3

UMAPPlot(Chrimsonposcells)
ggsave("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_cluster_umap.png", height = 4, width = 5)
UMAPPlot(Chrimsonposcells, group.by='seurat_clusters_two_groups')
ggsave("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_metacluster_umap.png", height = 4, width = 4)

## DE between these two groups
DE_Chrimsonposcells_two_groups <- give_top_logf_genes(Seurat::FoldChange(Chrimsonposcells, ident.1='TRUE', ident.2='FALSE', group.by='seurat_clusters_two_groups'))
DE_Chrimsonposcells_two_groups
xtable::xtable(DE_Chrimsonposcells_two_groups)
FeaturePlot(Chrimsonposcells, features = rownames(DE_Chrimsonposcells_two_groups))
ggsave("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_metaclustrer_featureplot.png", height = 6, width = 8)

## DE between light on and light off only in Chrimson-positive cells
Chrimsonposcells@meta.data$stim_light <- (Chrimsonposcells@meta.data$stim %in% c('G1', 'G3'))
Chrimsonposcells@meta.data$stim_light[Chrimsonposcells@meta.data$stim_light] <- 'Lights on'
Chrimsonposcells@meta.data$stim_light[Chrimsonposcells@meta.data$stim_light == 'FALSE'] <- 'Lights off'
table(Chrimsonposcells@meta.data$stim_light)
DE_Chrimsonposcells_lights <- Seurat::FoldChange(Chrimsonposcells, ident.1='Lights on', ident.2='Lights off', group.by='stim_light')
DE_Chrimsonposcells_lights
DE_Chrimsonposcells_lights_topgenes <- give_top_logf_genes(DE_Chrimsonposcells_lights)
xtable::xtable(DE_Chrimsonposcells_lights_topgenes)
FeaturePlot(Chrimsonposcells, features = rownames(DE_Chrimsonposcells_lights_topgenes))
ggsave("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_lightONvsOFF_featureplot.png", height = 6, width = 8)

ggplot(DE_Chrimsonposcells_lights, aes(x=avg_log2FC, y=pct.1))+
  geom_point(alpha=0.01)+lims(x=c(-0.25, 0.25))
EnhancedVolcano::EnhancedVolcano(DE_Chrimsonposcells_lights, x='avg_log2FC', y='pct.1', lab=rownames(DE_Chrimsonposcells_lights))

## DE of specific clusters

DE_Chrimsonposcells_lights_per_cluster <- give_cluster_specific_DE(Chrimsonposcells,
                                                                   cluster_name='seurat_clusters',
                                                                   ident.1='Lights on', ident.2='Lights off', group.by='stim_light')
DE_Chrimsonposcells_lights_per_cluster_topgenes <- lapply(DE_Chrimsonposcells_lights_per_cluster, give_top_logf_genes, include_separation_row=F, n=4)
DE_Chrimsonposcells_lights_per_cluster_topgenes <- lapply(DE_Chrimsonposcells_lights_per_cluster_topgenes, function(i) data.frame(gene=rownames(i), avg_log2FC=i[,'avg_log2FC']))
DE_Chrimsonposcells_lights_per_cluster_topgenes <- melt(DE_Chrimsonposcells_lights_per_cluster_topgenes)
table(DE_Chrimsonposcells_lights_per_cluster_topgenes$variable)
DE_Chrimsonposcells_lights_per_cluster_topgenesdcast <- dcast(DE_Chrimsonposcells_lights_per_cluster_topgenes, L1~gene, value.var = "value")
DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[is.na(DE_Chrimsonposcells_lights_per_cluster_topgenesdcast)] <- 0 ## semi-controversial

bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)

DE_Chrimsonposcells_lights_per_cluster_topgenes$gene <- factor(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene, levels=names(sort(table(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene))))
ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes, aes(y=gene, x=L1, fill=value))+
  scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98")+
  geom_tile()+ggtitle('Differentially expressed genes between lights ON/OFF per cell cluster')+
  theme_bw()+labs(x='Cluster', y='Gene')
ggsave("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_lightONvsOFF_percluster_DEgenes.png", height = 11, width = 8)

## to find clustering
select_two_or_more_active <- function(i){
  i[,colSums(i>0) > 2]
}
pdf("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_lightONvsOFF_percluster_DEgenes_subsetgenes_cluster.pdf", height = 3, width = 4)
print(pheatmap::pheatmap(t(as(select_two_or_more_active(DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1]), 'matrix'))))
dev.off()

install.packages('ggjoy')
library(ggjoy)
table(DE_Chrimsonposcells_lights_per_cluster_topgenes$variable)
ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes %>% dplyr::filter(gene %in% c('sr', 'Hr38', 'CG14186', 'cbt', 'CG46385')), aes(y=gene, x=value))+geom_joy() + theme_joy() 
ggsave("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_lightONvsOFF_percluster_DEgenes_subsetgenes_logfoldchange.pdf", height = 3, width = 4)

## summary of most important genes across clusters
ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes[DE_Chrimsonposcells_lights_per_cluster_topgenes$gene %in% tail(levels(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene)),], aes(y=gene, x=L1, fill=value))+
  scale_fill_gradientn(colours=c(bl,"white", re))+
  geom_tile()+ggtitle('Differentially expressed genes between lights ON/OFF per cell cluster\n(most shared genes)')+
  theme_bw()+labs(x='Cluster', y='Gene')
ggsave("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_lightONvsOFF_percluster_DEgenes_subsetgenes.png", height = 4, width = 8)

splitUMAPPlot(Chrimsonposcells, group.by='stim')
splitUMAPPlot(Chrimsonposcells, group.by='stim_light')

## UMAP of lights on/off coloured by feature
splitUMAPPlot(Chrimsonposcells, group.by='stim_light', colour_by_feature='sr')
ggsave("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_lightONvsOFF_sr.png", height = 2.2, width = 4)
splitUMAPPlot(Chrimsonposcells, group.by='stim_light', colour_by_feature='Hr38')
ggsave("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_lightONvsOFF_Hr38.png", height = 2.2, width = 4)
splitUMAPPlot(Chrimsonposcells, group.by='stim_light', colour_by_feature='Hr38', dim_red = "tsne")
ggsave("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_lightONvsOFF_Hr38_TSNE.png", height = 2.2, width = 4)
splitUMAPPlot(Chrimsonposcells, group.by='stim_light', colour_by_feature='sr', dim_red = "tsne")
ggsave("~/projects/lighting/3_results/chrimsonpos/chrimsonpos_lightONvsOFF_sr_TSNE.png", height = 2.2, width = 4)

## percentages of genes which are DE, in each cluster, including all cells, or separating by Chrimson cells

## --------------------------------------------------------- ##
## Analyses without the transposable elements
# load(file = "~/projects/lighting/data/robjects/image_DEA_Charly.RData")

FEB2019_noTE <- list(FEB2019_G1_rep1=FEB2019_G1_rep1,
                     FEB2019_G1_rep2=FEB2019_G1_rep2,
                     FEB2019_G2_rep1=FEB2019_G2_rep1,
                     FEB2019_G2_rep2=FEB2019_G2_rep2,
                     FEB2019_G3_rep1=FEB2019_G3_rep1,
                     FEB2019_G3_rep2=FEB2019_G3_rep2,
                     FEB2019_G4_rep1=FEB2019_G4_rep1,
                     FEB2019_G4_rep2=FEB2019_G4_rep2)
## remove TE
subsets_noTE <- lapply(FEB2019_noTE, function(i) which(!grepl("^TE-", rownames(i@assays$RNA@counts))))
FEB2019_noTE <- lapply(names(FEB2019_noTE), function(i) subset(FEB2019_noTE[[i]], features = subsets_noTE[[i]]))

## normalise data
FEB2019_noTE <- lapply(FEB2019_noTE, NormalizeData)

## find variable features
FEB2019_noTE <- lapply(FEB2019_noTE, FindVariableFeatures)
FEB2019_noTE <- FindIntegrationAnchors(object.list = FEB2019_noTE, dims = 1:60)
FEB2019_noTE <- IntegrateData(anchorset = FEB2019_noTE, dims = 1:60)
DefaultAssay(FEB2019_noTE) <- "integrated"
FEB2019_noTE <- add_dimred(FEB2019_noTE)

saveRDS(FEB2019_noTE, "~/projects/lighting/data/robjects/FEB2019_noTE.RDS")
# FEB2019_noTE <- readRDS("~/projects/lighting/data/robjects/FEB2019_noTE.RDS")

FEB2019_noTE <- FindClusters(FEB2019_noTE, resolution = 0.1)
FEB2019_noTE <- FindClusters(FEB2019_noTE, resolution = 0.01)
FEB2019_noTE <- FindClusters(FEB2019_noTE, resolution = 0.001)

## find latest resolution
FEB2019_noTE@meta.data$integrated_snn_res.4
FEB2019_noTE@meta.data$integrated_snn_res.0.1

FEB2019_noTE_markers_0p001 <- FindAllMarkers(FEB2019_noTE, logfc.threshold = 0.8, assay = "RNA")

FEB2019_noTE_markers_0p001
# Seurat::DoHeatmap(FEB2019_noTE) ## it takes a long time ## prohibitive
# Seurat::DoHeatmap(FEB2019_noTE, features = unique(FEB2019_noTE_markers_0p001$gene)) ## prohibitive
VlnPlot(FEB2019_noTE, features = unique(FEB2019_noTE_markers_0p001$gene))
# pheatmap::pheatmap(as(FEB2019_noTE@assays$RNA[unique(FEB2019_noTE_markers_0p001$gene),], 'matrix')) ## prohibitive


pheatmap::pheatmap(FEB2019_noTE_markers_0p001[,c('pct.1', 'pct.2')])

FeaturePlot(FEB2019_noTE, features='Chrimson')
FeaturePlot(FEB2019_noTE, features='lncRNA:roX1')
FeaturePlot(FEB2019_noTE, features='Chrimson', reduction = "tsne")
UMAPPlot(FEB2019_noTE, group.by='integrated_snn_res.0.1')
UMAPPlot(FEB2019_noTE, group.by='integrated_snn_res.0.01') ## nice
TSNEPlot(FEB2019_noTE, group.by='integrated_snn_res.0.01') ## 
UMAPPlot(FEB2019_noTE, group.by='integrated_snn_res.0.001') ## bit weird
TSNEPlot(FEB2019_noTE, group.by='integrated_snn_res.0.001') ##
FeaturePlot(FEB2019_noTE, features='Chrimson', reduction = "tsne")
UMAPPlot(FEB2019_noTE)
