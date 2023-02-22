## analyse DEG from Charly
## using some of his code

rm(list = ls())
setwd("/t1-data/project/cncb/shared/proj002/analyses/2020/10kCells/")

library(Seurat)
library(ggplot2)
library(cowplot)
library(plotly)
library(reshape2)
library(ggjoy)
library(jcolors)
library(ggrepel)
library(SingleR)
theme_set(theme_cowplot())

source("~/projects/lighting/2_analyses/helper_functions.R")

##---------------------------------------------------------------------------------#
## Biomarkers
biomarkers_list <- readRDS("/t1-data/project/sims-lab/lmorrill/robjects/lightning/biomarkers_list.RDS")

##---------------------------------------------------------------------------------#
## Read in data

# input_objs <- 'CharlyCellRanger'
input_objs <- 'LenaCellRanger'

if(input_objs=='CharlyCellRanger'){
  
  Read10X("./FEB2019_G1_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G1_rep1.data
  Read10X("./FEB2019_G1_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G1_rep2.data
  Read10X("./FEB2019_G2_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G2_rep1.data
  Read10X("./FEB2019_G2_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G2_rep2.data
  Read10X("./FEB2019_G3_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G3_rep1.data
  Read10X("./FEB2019_G3_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G3_rep2.data
  Read10X("./FEB2019_G4_rep1/outs/filtered_feature_bc_matrix") -> FEB2019_G4_rep1.data
  Read10X("./FEB2019_G4_rep2/outs/filtered_feature_bc_matrix") -> FEB2019_G4_rep2.data
}else if(input_objs == 'LenaCellRanger'){
  folderLenaCellRanger <- '/project/sims-lab/lmorrill/data/lightning/alignment/cellrangerOUT/dmel649CORRECTED/'
  Read10X(paste0(folderLenaCellRanger, "G1_rep1/G1_rep1/outs/filtered_feature_bc_matrix/")) -> FEB2019_G1_rep1.data
  Read10X(paste0(folderLenaCellRanger, "G1_rep2/G1_rep2/outs/filtered_feature_bc_matrix/")) -> FEB2019_G1_rep2.data
  Read10X(paste0(folderLenaCellRanger, "G2_rep1/G2_rep1/outs/filtered_feature_bc_matrix/")) -> FEB2019_G2_rep1.data
  Read10X(paste0(folderLenaCellRanger, "G2_rep2/G2_rep2/outs/filtered_feature_bc_matrix/")) -> FEB2019_G2_rep2.data
  Read10X(paste0(folderLenaCellRanger, "G3_rep1/G3_rep1/outs/filtered_feature_bc_matrix/")) -> FEB2019_G3_rep1.data
  Read10X(paste0(folderLenaCellRanger, "G3_rep2/G3_rep2/outs/filtered_feature_bc_matrix/")) -> FEB2019_G3_rep2.data
  Read10X(paste0(folderLenaCellRanger, "G4_rep1/G4_rep1/outs/filtered_feature_bc_matrix/")) -> FEB2019_G4_rep1.data
  Read10X(paste0(folderLenaCellRanger, "G4_rep2/G4_rep2/outs/filtered_feature_bc_matrix/")) -> FEB2019_G4_rep2.data
  
}else{
  stop()
}

FEB2019_G1_rep1 <- CreateSeuratObject(FEB2019_G1_rep1.data, project = "FEB2019_G1_rep1", min.cells = 3, min.features = 200)
FEB2019_G1_rep2 <- CreateSeuratObject(FEB2019_G1_rep2.data, project = "FEB2019_G1_rep2", min.cells = 3, min.features = 200)
FEB2019_G2_rep1 <- CreateSeuratObject(FEB2019_G2_rep1.data, project = "FEB2019_G2_rep1", min.cells = 3, min.features = 200)
FEB2019_G2_rep2 <- CreateSeuratObject(FEB2019_G2_rep2.data, project = "FEB2019_G2_rep2", min.cells = 3, min.features = 200)
FEB2019_G3_rep1 <- CreateSeuratObject(FEB2019_G3_rep1.data, project = "FEB2019_G3_rep1", min.cells = 3, min.features = 200)
FEB2019_G3_rep2 <- CreateSeuratObject(FEB2019_G3_rep2.data, project = "FEB2019_G3_rep2", min.cells = 3, min.features = 200)
FEB2019_G4_rep1 <- CreateSeuratObject(FEB2019_G4_rep1.data, project = "FEB2019_G4_rep1", min.cells = 3, min.features = 200)
FEB2019_G4_rep2 <- CreateSeuratObject(FEB2019_G4_rep2.data, project = "FEB2019_G4_rep2", min.cells = 3, min.features = 200)

rm(list=ls(pattern=".data"))

all_samples <- list(FEB2019_G1_rep1,
                    FEB2019_G1_rep2,
                    FEB2019_G2_rep1,
                    FEB2019_G2_rep2,
                    FEB2019_G3_rep1,
                    FEB2019_G3_rep2,
                    FEB2019_G4_rep1,
                    FEB2019_G4_rep2)

sample_names <- apply(expand.grid(paste0('G', 1:4), paste0('_rep', 1:2)), 1, paste0, collapse='')
sample_names
names(all_samples) <- sample_names

rm(FEB2019_G1_rep1)
rm(FEB2019_G1_rep2)
rm(FEB2019_G2_rep1)
rm(FEB2019_G2_rep2)
rm(FEB2019_G3_rep1)
rm(FEB2019_G3_rep2)
rm(FEB2019_G4_rep1)
rm(FEB2019_G4_rep2)

## Change name of features, if necessary
if(input_objs == 'LenaCellRanger'){
  
  name_conversion_file <- give_name_conversion_file()
  for(i in 1:length(all_samples)){
    stopifnot(rownames(all_samples[[i]]@assays$RNA@data) == all_samples[[i]]@assays$RNA@counts@Dimnames[[1]])
    gene_names_replace <- name_conversion_file[,2][match(all_samples[[i]]@assays$RNA@counts@Dimnames[[1]], name_conversion_file[,1])]
    table(is.na(gene_names_replace))
    stopifnot(is.na(gene_names_replace) == 0)
    all_samples[[i]]@assays$RNA@counts@Dimnames[[1]] <- gene_names_replace ## first change
    rownames(all_samples[[i]]@assays$RNA@data) <- gene_names_replace ## second change
  }
  
  # name_conversion_file[,2][match('FBgn0024733', name_conversion_file[,1])]
  # head(cbind.data.frame(gene_names_replace, FB=FEB2019_G1_rep1@assays$RNA@counts@Dimnames[[1]]))
}

FeatureScatters <- list()
FeatureScatters_2 <- list()
for(i in 1:length(all_samples)){
  all_samples[[i]][["percent.mt"]] <- PercentageFeatureSet(all_samples[[i]], pattern = "^mt:")
  FeatureScatters[[i]] <- FeatureScatter(all_samples[[i]], feature1 = "nCount_RNA", feature2 = "percent.mt")
  FeatureScatters_2[[i]] <- FeatureScatter(all_samples[[i]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
}

for(i in 1:length(all_samples)){
  print(sum(all_samples[[i]][["percent.mt"]]))
}

CombinePlots(FeatureScatters)
CombinePlots(FeatureScatters_2)

for(i in 1:length(all_samples)){
  all_samples[[i]] <- subset(all_samples[[i]], subset = nCount_RNA < 20000 & percent.mt < 15)
}

## note we are not subsetting cells with very few reads - possibly, those are subsetted by <cellranger count> already

# plot(density(FEB2019_G1_rep1$nCount_RNA))
# plot(density(log(FEB2019_G1_rep1$nCount_RNA))) ## bimodality seems to indicate a population of cells with very few reads
# plot(density(log(FEB2019_G1_rep1$nFeature_RNA)))
# plot((log(FEB2019_G1_rep1$nFeature_RNA)), log(FEB2019_G1_rep1$nCount_RNA), col=alpha('black', 0.02))

for(i in 1:length(all_samples)){
  all_samples[[i]][["RNA"]]@meta.features <- data.frame(row.names = rownames(all_samples[[i]][["RNA"]]))
  all_samples[[i]] <- NormalizeData(all_samples[[i]])
  all_samples[[i]] <- FindVariableFeatures(all_samples[[i]]) ## https://github.com/satijalab/seurat/issues/2317
}

anchors <- FindIntegrationAnchors(object.list = all_samples, dims = 1:60)
combined_dataset <- IntegrateData(anchorset = anchors, dims = 1:60)

DefaultAssay(combined_dataset) <- "integrated"
combined_dataset <- ScaleData(combined_dataset)
combined_dataset <- RunPCA(combined_dataset, npcs = 60)
combined_dataset <- RunUMAP(combined_dataset, reduction = "pca", dims = 1:60)
combined_dataset <- FindNeighbors(combined_dataset, reduction = "pca", dims = 1:60)
combined_dataset <- FindClusters(combined_dataset, resolution = 4)
combined_dataset <- RunTSNE(combined_dataset, reduction = "pca", dims = 1:60)

if(input_objs == "LenaCellRanger"){
  combined_dataset$stim <- combined_dataset$orig.ident
}

combined_dataset$stim <- gsub("_.*", "", gsub('FEB2019_', '', combined_dataset$stim))
combined_dataset$stim

DefaultAssay(combined_dataset) <- "RNA"
MARKERS.combined_dataset <- FindAllMarkers(combined_dataset, logfc.threshold = 0.8, assay = "RNA")
# save(MARKERS.combined_dataset, file = "./FEB2019.MARKERS.allgroups.Robj")
system(paste0("mkdir -p ~/projects/lighting/data/robjects/", input_objs))
save(MARKERS.combined_dataset, file = paste0("~/projects/lighting/data/robjects/", input_objs, "/image_Markers.RData"))


ALL.average.expression <- AverageExpression(combined_dataset, verbose = F)

Idents(combined_dataset) <- combined_dataset$integrated_snn_res.4
# KCs.cell <- subset(combined_dataset, idents = c("3", "7", "29", "32", "36", "78", "90", "100", "112"))
# Idents(KCs.cell) <- "stim"
# avg.KCs.cell <- log1p(AverageExpression(KCs.cell, verbose = F)$RNA)
# avg.KCs.cell$gene <- rownames(avg.KCs.cell)
# p1 <- ggplot(avg.KCs.cell, aes(G3, G2, name = gene)) + geom_point() + ggtitle("G3 vs. G2")
# ggplotly(p1)
# p2 <- ggplot(avg.KCs.cell[telist,], aes(G3, G2, name = gene)) + geom_point() + ggtitle("G3 vs. G2")

VlnPlot(combined_dataset, "Dsk", group.by = "stim", split.by = "orig.ident", pt.size = 0.1, slot = "counts", assay = "RNA", y.max = 10)


# for (i in telist) {
#   jpeg(file = paste("~/Dropbox/CloudDesktop/TEplots/TE_count_", i, ".jpeg", sep=""))
#   plot1 <- VlnPlot(combined_dataset, idents = c(27), group.by = "orig.ident",features = as.character(i), slot = "counts", assay = "RNA")
#   plot(plot1)
#   dev.off()
# }

save.image(file = paste0("~/projects/lighting/data/robjects/", input_objs, "/image_DEA.RData"))

## --------------------------------------------------------- ##
## FROM NOW ON: Lena
## --------------------------------------------------------- ##
# load(file = "~/projects/lighting/data/robjects/image_DEA_Charly.RData")

## --------------------------------------------------------- ##
## Cluster annotation

## first, at very low resolution
table(combined_dataset$integrated_snn_res.4)
table(combined_dataset$stim)

DefaultAssay(combined_dataset) <- 'integrated'
combined_dataset <- FindClusters(combined_dataset, resolution = 0.1)

UMAPPlot(combined_dataset, group.by='integrated_snn_res.0.1')

DefaultAssay(combined_dataset) <- "RNA"
MARKERS_lower_granularity <- FindAllMarkers(combined_dataset, logfc.threshold = 0.8, assay = "RNA")

## select top gene for each cluster. the df from FindAllMarkers is already sorted by p-value and LFC
MARKERS_lower_granularity
MARKERS_lower_granularity_top1 <- MARKERS_lower_granularity[sapply(unique(MARKERS_lower_granularity$cluster), function(i) which(MARKERS_lower_granularity$cluster == i)[1]),]
MARKERS_lower_granularity_top1 <- MARKERS_lower_granularity_top1[,c('avg_log2FC', 'gene', 'cluster')]
MARKERS_lower_granularity_top1_dcast <- na_to_zeros(firstcol_to_rownames(dcast(MARKERS_lower_granularity_top1, gene~cluster, value.var = 'avg_log2FC')))

MARKERS_lower_granularity_top1[order(MARKERS_lower_granularity_top1$gene),]
pheatmap::pheatmap(MARKERS_lower_granularity_top1_dcast)

ßif(input_objs =='CharlyCellRanger'){
ggplot(MARKERS_lower_granularity_top1, aes(x=as.numeric(cluster), y = factor(gene, levels=rev(gtools::mixedsort(gene))), fill=avg_log2FC))+
  geom_tile()+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 1, label = "Glutamatergic neurons", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 2, label = "Vestigial+ cells", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 3, label = "Unknown", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 4, label = "Cholinergic neurons", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 5, label = "Possibly CRZ neurons", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 6, label = "Ellipsoid body", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 7, label = "Dorsal tract", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 8, label = "Unknown", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 9, label = "Histamine", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 10, label = "Unknown", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 11, label = "α'/β' mushroom body neurons", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 12, label = "Transmedullary neurons", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 13, label = "PNs", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 14, label = "Ventral PNs", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 15, label = "Unknown", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 16, label = "Sp1+ neural progeny", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 17, label = "Dorsal clock neurons", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 18, label = "Dopamine neurons", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 19, label = "Unknown", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 20, label = "Unknown", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 21, label = "Serotonergic neurons", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 22, label = "Unknown", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 23, label = "Unknown", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 24, label = "Unknown", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 25, label = "Fat body", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 26, label = "L4 and L5 lamina neurons and or Mi1 medulla neurons", hjust = 0)+
  annotate("text", x = length(MARKERS_lower_granularity_top1$cluster)+10, y = 27, label = "Astrocytes", hjust = 0)+
  xlim(c(0, length(MARKERS_lower_granularity_top1$cluster)+20))+
  theme_bw()


annotation_clusters_0p1 <- cbind.data.frame(cluster=MARKERS_lower_granularity_top1$cluster[order(MARKERS_lower_granularity_top1$gene,
                                                                                                 decreasing = T)],
                                            gene=sort(MARKERS_lower_granularity_top1$gene, decreasing = T),
  annotation=c("Glutamatergic neurons", "Vestigial+ cells", "Unknown", "Cholinergic neurons", "Possibly CRZ neurons", "Ellipsoid body",
    "Dorsal tract", "Unknown", "Histamine", "Unknown", "α'/β' mushroom body neurons", "Transmedullary neurons", "PNs",
    "Ventral PNs", "Unknown", "Sp1+ neural progeny", "Dorsal clock neurons", "Dopamine neurons", "Unknown", "Unknown",
    "Serotonergic neurons", "Unknown", "Unknown", "Unknown", "Fat body", "L4 and L5 lamina neurons and or Mi1 medulla neurons",
    "Astrocytes"))
annotation_clusters_0p1


print(xtable::xtable(cbind.data.frame(annotation=annotation_clusters_0p1$annotation[match(unique(MARKERS_lower_granularity_top1$cluster), 
                                                                                          annotation_clusters_0p1$cluster)],
MARKERS_lower_granularity_top1[,c(3,2,1)])), include.rownames=FALSE)

## table of multiple other markers
MARKERS_lower_granularity_top10 <- MARKERS_lower_granularity[sapply(unique(MARKERS_lower_granularity$cluster), function(i) which(MARKERS_lower_granularity$cluster == i)[1:10]),]
## into table

print(xtable::xtable(cbind.data.frame(cluster=unique(MARKERS_lower_granularity_top10$cluster),
                                      annotation=annotation_clusters_0p1$annotation[match(unique(MARKERS_lower_granularity_top10$cluster), annotation_clusters_0p1$cluster)],
                 markers=sapply(unique(MARKERS_lower_granularity_top10$cluster), function(i) paste0(sort(MARKERS_lower_granularity_top10$gene[MARKERS_lower_granularity_top10$cluster == i]), collapse = ', '))),
), include.rownames = F)
}

## at even lower resolution, what are the clusters? Glia vs other?
DefaultAssay(combined_dataset) <- 'integrated'
combined_dataset <- FindClusters(combined_dataset, resolution = 0.01)
DefaultAssay(combined_dataset) <- "RNA"
MARKERS_lowest_granularity <- FindAllMarkers(combined_dataset, logfc.threshold = 0.8, assay = "RNA")
UMAPPlot(combined_dataset, group.by='integrated_snn_res.0.01')

table(combined_dataset$integrated_snn_res.0.1,
combined_dataset$integrated_snn_res.0.01)

pheatmap::pheatmap(normalise_rw(table(combined_dataset$integrated_snn_res.0.1,
                         combined_dataset$integrated_snn_res.0.01)))

## cluster annotation using SingleR
# library(celldex)
# ref <- BlueprintEncodeData()
# 
# pred <- SingleR(test=combined_dataset, ref=ref, labels=ref$label.main)

## add top marker to each cluster in UMAP -- facet by cluster

umap_facets_with_topmarker <- function(seurat_obj, markers, seurat_name_clusters){
  reduction_df <- data.frame(Reductions(seurat_obj, 'umap')@cell.embeddings)
  reduction_df$clusters <- seurat_obj[[seurat_name_clusters]][,1]
  reduction_df$topmarker <- NA
  reduction_df[sapply(unique(reduction_df$clusters), function(i) which(reduction_df$cluster == i)[1]), 'topmarker'] <- markers$gene[sapply(unique(reduction_df$clusters), function(i) which(markers$cluster == i)[1])]
  table(reduction_df$topmarker)
  
  ggplot(reduction_df, aes(x=UMAP_1, y=UMAP_2, label=topmarker))+geom_point()+facet_wrap(.~clusters)+theme_bw()+
    geom_label(aes(x = 0, y=0))
  
}

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

# col_vector <- sapply(1:100, function(i) paste0(c('#', sample(c(0:9, LETTERS), size = 6)), collapse = ''))
col_vector <- sample(colours(), size = 100)
# col_vector <- readRDS("~/small_practical_robjects/col_vector1.RDS")
# saveRDS(col_vector, "~/small_practical_robjects/col_vector3.RDS")

umap_single_facet_with_topmarker <- function(seurat_obj, markers, seurat_name_clusters,density_map=F){
  reduction_df <- data.frame(Reductions(seurat_obj, 'umap')@cell.embeddings)
  reduction_df$clusters <- seurat_obj[[seurat_name_clusters]][,1]
  reduction_df$topmarker <- NA
  reduction_df <- rbind(reduction_df,
  cbind.data.frame(UMAP_1 = sapply(unique(reduction_df$clusters), function(i) mean(reduction_df$UMAP_1[which(reduction_df$cluster == i)])),
                   UMAP_2= sapply(unique(reduction_df$clusters), function(i) mean(reduction_df$UMAP_2[which(reduction_df$cluster == i)])),
                   clusters = unique(reduction_df$clusters),
                   topmarker = markers$gene[sapply(unique(reduction_df$clusters), function(i) which(markers$cluster == i)[1])])
  )
  
  if(density_map){
    ggplot(reduction_df, aes(x=UMAP_1, y=UMAP_2, label=topmarker, col=clusters))+geom_density_2d()+theme_bw()+
      geom_label_repel()+
      scale_color_manual(values = col_vector)    
  }else{
    ## points
    ggplot(reduction_df, aes(x=UMAP_1, y=UMAP_2, label=topmarker, col=clusters))+geom_point()+theme_bw()+
      geom_label_repel()+
      scale_color_manual(values = col_vector)
  }
  
}

umap_single_facet_with_topmarker(seurat_obj = combined_dataset, markers = MARKERS_lowest_granularity, seurat_name_clusters = 'integrated_snn_res.0.01')
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/umap_clusters_lowest_granularity_markers.pdf"), height = 6, width = 6)

umap_facets_with_topmarker(seurat_obj = combined_dataset, markers = MARKERS_lowest_granularity, seurat_name_clusters = 'integrated_snn_res.0.01')
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/umap_clusters_lowest_granularity_markers_facets.pdf"), height = 8, width = 8)

umap_single_facet_with_topmarker(seurat_obj = combined_dataset, markers = MARKERS_lower_granularity, seurat_name_clusters = 'integrated_snn_res.0.1')
umap_single_facet_with_topmarker(seurat_obj = combined_dataset, markers = MARKERS_lower_granularity, seurat_name_clusters = 'integrated_snn_res.0.1', density_map=T)

## too many clusters
# umap_single_facet_with_topmarker(seurat_obj = combined_dataset, markers = MARKERS.combined_dataset, seurat_name_clusters = 'integrated_snn_res.4')

# seurat_obj = combined_dataset
# markers = MARKERS_lowest_granularity
# seurat_name_clusters = 'integrated_snn_res.0.01'

## go term of markers of each cluster
MARKERS_lowest_granularity[MARKERS_lowest_granularity$cluster == 0,'gene']

## --------------------------------------------------------- ##


## --------------------------------------------------------- ##
## what are the Chrimson positive cells?
#' es, absolutely: It’s mainly a bunch of positively reinforcing dopamine neurons, plus a group of off-target 
#' cholinergic neurons (plus cells in the eye, but they are not part of the data set) . The driver we used is 0273-Gal4

combined_dataset['Chrimson']
plotdensity <- function(arg_vec, log2p1lusimp=T, value_imput=1){
  x <- (reshape2::melt(arg_vec))
  if(log2p1lusimp){
    ggplot(x, aes(x=value+value_imput))+geom_density()+scale_x_continuous(trans = "log10")
  }else{
    ggplot(x, aes(x=value))+geom_density()+scale_x_continuous(trans = "log2")
  }
}
plot(density(as.vector(combined_dataset@assays$RNA['Chrimson'])))
plotdensity(as.vector(combined_dataset@assays$RNA['Chrimson']), value_imput = 0.0000001)
plotdensity(as.vector(combined_dataset@assays$RNA['Chrimson']), value_imput = 0)
combined_dataset@assays$integrated['Chrimson']
max(as.vector(combined_dataset@assays$RNA['Chrimson']))
plotdensity(as.vector(combined_dataset@assays$RNA['Chrimson']), value_imput = 0)


## --------------------------------------------------------- ##
## umaps
system(paste0("mkdir -p ~/projects/lighting/3_results/", input_objs))

# save(MARKERS.combined_datasetfile = paste0("~/projects/lighting/3_results/", input_objs, "))

UMAPPlot(combined_dataset, group.by='stim')
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/stim_umap.png"), height = 4, width = 4)

splitUMAPPlot(combined_dataset, group.by='stim')
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/stim_umap_split.png"), height = 6, width = 8)

FeaturePlot(combined_dataset, features = "Chrimson")
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/chrimsom_umap.png"), height = 4, width = 4)

FeaturePlot(combined_dataset, features = "Chrimson", reduction = "tsne")
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/chrimsom_tsne.png"), height = 4, width = 4)

DimPlot(combined_dataset, reduction = "umap")+ NoLegend()
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/clusters_umap.pn"), height = 4, width = 4)

DimPlot(combined_dataset, reduction = "tsne")+ NoLegend()
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/clusters_tsne.png"), height = 4, width = 4)

combined_dataset@assays

FeaturePlot(combined_dataset, features = biomarkers_list$`gene name`)

VlnPlot(combined_dataset, features  = biomarkers_list$`gene name`)

# aggregate(x = as(combined_dataset@assays$RNA@counts['Vmat',], 'vector'), by = lapply(unique( combined_dataset$seurat_clusters), function(i)  combined_dataset$seurat_clusters == i ), FUN='mean')
Vmat_across_clusters <- aggregate(x = as(combined_dataset@assays$RNA@counts['Vmat',], 'vector'), by = list(clusters=combined_dataset$seurat_clusters), FUN='mean')
ggplot(Vmat_across_clusters, aes(x=factor(clusters, levels=clusters[order(x)]), y=x))+geom_point()+geom_line()

biomarkers_list_vec <- c(biomarkers_list$`gene name`, 'Vmat', 'DAT', 'Frq1', 'otp', 'Lmx1a', 'Tdc2', 'Tbh', 'Tbh', 'kek1', 'hth', 'mirr', 'vvl', 'ey', 'Lim3')
## aggregate by mean
# biomarkers_across_clusters <- sapply(biomarkers_list_vec, function(biomarker_it){
#   aggregate(x = as(combined_dataset@assays$RNA@counts[biomarker_it,], 'vector'), by = list(clusters=combined_dataset$seurat_clusters), FUN='mean')[,2]
# })
## aggregate by median
## Charlys and my dataset have subtle differences in the gene names
biomarkers_list_vec <- unique(c(biomarkers_list_vec, biomarkers_list$`other name`[!(biomarkers_list_vec %in% rownames(combined_dataset@assays$RNA@counts))]))
biomarkers_list_vec <- biomarkers_list_vec[(biomarkers_list_vec %in% rownames(combined_dataset@assays$RNA@counts))]

biomarkers_across_clusters <- sapply(biomarkers_list_vec, function(biomarker_it){
  aggregate(x = as(combined_dataset@assays$RNA@counts[biomarker_it,], 'vector'), by = list(clusters=combined_dataset$seurat_clusters), FUN='median')[,2]
})
pheatmap::pheatmap(biomarkers_across_clusters)

pheatmap::pheatmap(remove_cols_with_na(normalise_cl(biomarkers_across_clusters)))

FeaturePlot(combined_dataset, features=biomarkers_list$`gene name`[grepl('KC', biomarkers_list$marker)]) ## KC markers

pheatmap::pheatmap(ALL.average.expression$RNA[biomarkers_list$`gene name`[grepl('KC', biomarkers_list$marker)],]) ## KC markers
pheatmap::pheatmap(scale(ALL.average.expression$RNA[biomarkers_list$`gene name`[grepl('KC', biomarkers_list$marker)],])) ## KC markers. scaled data. no clear clustering of KC vs non-KC?
pheatmap::pheatmap(scale(center = T, ALL.average.expression$RNA[biomarkers_list_vec,]))
pheatmap::pheatmap(scale(center = T, ALL.average.expression$integrated[remove_na(match(biomarkers_list_vec, rownames(ALL.average.expression$integrated))),]))

dim(ALL.average.expression$RNA)
dim(ALL.average.expression$integrated)

# Seurat::clust (combined_dataset)

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
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/chrimsom_umap_", i, ".png"), height = 4, width = 4)

}

rownames(combined_dataset$nCount_RNA)

Seurat::Assays(combined_dataset)

combined_dataset@assays$RNA@counts['Chrimson',]
plot(density(combined_dataset@assays$RNA@counts['Chrimson',]))
plot(density(log(combined_dataset@assays$RNA@counts['Chrimson',])))
plot(density(log(1+combined_dataset@assays$RNA@counts['Chrimson',])))
plot(density(as.vector(FetchData(combined_dataset, 'Chrimson'))[[1]]))
table(as.vector(FetchData(combined_dataset, 'Chrimson'))[[1]] > 0)

## hard threshold of at least one count in Chrimson (i.e. probably many false negatives)
table(combined_dataset@assays$RNA@counts['Chrimson',] > 0)

##---------------------------------------------------------------------------------#
## Split by Chrimson

Chrimsonpos <- combined_dataset@assays$RNA@counts['Chrimson',] > 0
table(Chrimsonpos)

Chrimsonposcells <- subset(combined_dataset, cells = which(Chrimsonpos))
Chrimsonnegcells <- subset(combined_dataset, cells = which(!Chrimsonpos))
Chrimsonposcells ## 1899 samples
Chrimsonposcells <- FindVariableFeatures(Chrimsonposcells)
Chrimsonposcells <- add_dimred(Chrimsonposcells)
Chrimsonnegcells <- FindVariableFeatures(Chrimsonnegcells)
Chrimsonnegcells <- add_dimred(Chrimsonnegcells)

saveRDS(Chrimsonposcells, paste0("~/projects/lighting/data/robjects/", input_objs, "/Chrimsonposcells.RDS"))
saveRDS(Chrimsonnegcells, paste0("~/projects/lighting/data/robjects/", input_objs, "/Chrimsonnegcells.RDS"))

# Chrimsonposcells <- readRDS("~/projects/lighting/data/robjects/Chrimsonposcells.RDS")
# Chrimsonnegcells <- readRDS("~/projects/lighting/data/robjects/Chrimsonnegcells.RDS")

##---------------------------------------------------------------------------------#
## Biomarkers
## top markers for each of the clusters
xtable::xtable(MARKERS.combined_dataset[sapply(unique(MARKERS.combined_dataset$cluster), function(i) which(MARKERS.combined_dataset$cluster == i)[1]),])

UMAPPlot(Chrimsonposcells, group.by='stim')
UMAPPlot(Chrimsonposcells)
FeaturePlot(Chrimsonposcells, features='Chrimson')
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/chrimsonpos/chrimsonpos_UMAP_chrimson.pdf"), height = 3.5, width = 4)
FeaturePlot(Chrimsonposcells, features='Chrimson')
UMAPPlot(Chrimsonnegcells, group.by='stim')
UMAPPlot(Chrimsonnegcells)
TSNEPlot(Chrimsonposcells)
FeaturePlot(Chrimsonposcells, features='nCount_RNA')
FeaturePlot(Chrimsonposcells, features='nFeature_RNA')

table(Chrimsonposcells$stim)
table(Chrimsonposcells$stim_light)

plot(Chrimsonposcells$nFeature_RNA, Chrimsonposcells$nCount_RNA)

my_VlnPlot <- function(seurat_object, features, aggregation_fun='median', arg_ncol=4, logtrans=T){
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


## for each of the chrimson+ clusters, does any have move Chrimson than the others?
# my_VlnPlot(Chrimsonposcells, features = 'Chrimson') ## too slow? it doesn't work?
VlnPlot(Chrimsonposcells, features = 'Chrimson')

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
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/markers/",biomarkers_list[i,'gene name'], ".png"), height = 4, width = 4)
  FeaturePlot(Chrimsonposcells, features=biomarkers_list[i,'gene name'], reduction = "tsne")+
    labs(title = paste0(biomarkers_list[i,'gene name'],
                                                                                             "\n",
                                                                                             biomarkers_list[i,'marker']))
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/markers/TSNE_",biomarkers_list[i,'gene name'], ".png"), height = 4, width = 4)
  
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
ggsave("~/projects/lighting/3_results/", input_objs, "/chrimsonpos/chrimsonpos_cluster_umap.png", height = 4, width = 5)
UMAPPlot(Chrimsonposcells, group.by='seurat_clusters_two_groups')
ggsave("~/projects/lighting/3_results/", input_objs, "/chrimsonpos/chrimsonpos_metacluster_umap.png", height = 4, width = 4)

## DE between these two groups
DE_Chrimsonposcells_two_groups <- give_top_logf_genes(Seurat::FoldChange(Chrimsonposcells, ident.1='TRUE', ident.2='FALSE', group.by='seurat_clusters_two_groups'))
DE_Chrimsonposcells_two_groups
xtable::xtable(DE_Chrimsonposcells_two_groups)
FeaturePlot(Chrimsonposcells, features = rownames(DE_Chrimsonposcells_two_groups))
ggsave("~/projects/lighting/3_results/", input_objs, "/chrimsonpos/chrimsonpos_metaclustrer_featureplot.png", height = 6, width = 8)

## DE between light on and light off only in Chrimson-positive cells

give_DE_analysis_lightON_lightOFF <- function(dataset, dataset_name){
  dataset@meta.data$stim_light <- (dataset@meta.data$stim %in% c('G1', 'G3'))
  dataset@meta.data$stim_light[dataset@meta.data$stim_light] <- 'Lights on'
  dataset@meta.data$stim_light[dataset@meta.data$stim_light == 'FALSE'] <- 'Lights off'
  DE_Chrimsonposcells_lights <- Seurat::FoldChange(dataset, ident.1='Lights on', ident.2='Lights off', group.by='stim_light')
  DE_Chrimsonposcells_lights
  DE_Chrimsonposcells_lights_topgenes <- give_top_logf_genes(DE_Chrimsonposcells_lights)
  xtable::xtable(DE_Chrimsonposcells_lights_topgenes)
  FeaturePlot(dataset, features = rownames(DE_Chrimsonposcells_lights_topgenes))
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_featureplot.png"), height = 6, width = 8)

  ## volcano plot  
  # ggplot(DE_Chrimsonposcells_lights, aes(x=avg_log2FC, y=pct.1))+
  #   geom_point(alpha=0.01)+lims(x=c(-0.25, 0.25))
  # EnhancedVolcano::EnhancedVolcano(DE_Chrimsonposcells_lights, x='avg_log2FC', y='pct.1', lab=rownames(DE_Chrimsonposcells_lights))
  
  ## DE of specific clusters between conditions
  ## https://satijalab.org/seurat/archive/v3.0/de_vignette.html
  
  # FindMarkers(Chrimsonposcells, ident.1 = "Lights on", ident.2 = "Lights off", group.by = "stim_light")
  ## we are performing DE between conditions
  
  DE_Chrimsonposcells_lights_per_cluster <- give_cluster_specific_DE(dataset,
                                                                     cluster_name='seurat_clusters',
                                                                     ident.1='Lights on', ident.2='Lights off', group.by='stim_light')
  saveRDS(DE_Chrimsonposcells_lights_per_cluster, paste0("~/projects/lighting/data/robjects/", input_objs, "/", dataset_name, "_DEgenes_per_cluster.RDS"))
  
  # DE_Chrimsonposcells_lights_per_cluster_topgenes <- lapply(DE_Chrimsonposcells_lights_per_cluster, give_top_logf_genes, include_separation_row=F, n=4)
  DE_Chrimsonposcells_lights_per_cluster_topgenes <- lapply(DE_Chrimsonposcells_lights_per_cluster,function(i) i[i$p_val_adj <= 0.05,])
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
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes.png"), height = 11, width = 8)
  
  ## to find clustering
  select_two_or_more_active <- function(i){
    i[,colSums(i>0) > 2]
  }
  
  try(dev.off())
  pdf(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes_subsetgenes_cluster.pdf"), height = 3, width = 4)
  # print(pheatmap::pheatmap(t(as(select_two_or_more_active(DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1]), 'matrix'))))
  print(pheatmap::pheatmap(t(as((DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1]), 'matrix'))))
  dev.off()
  
  table(DE_Chrimsonposcells_lights_per_cluster_topgenes$variable)
  # ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes %>% dplyr::filter(gene %in% c('sr', 'Hr38', 'CG14186', 'cbt', 'CG46385')), aes(y=gene, x=value))+geom_joy() + theme_joy() 
  ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes, aes(y=gene, x=value))+geom_joy() + theme_joy() 
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes_subsetgenes_logfoldchange.pdf"), height = 3, width = 4)
  ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes[DE_Chrimsonposcells_lights_per_cluster_topgenes$gene %in% names(which(table(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene) > 2)),],
         aes(y=gene, x=value))+geom_joy() + theme_joy() 
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes_subsetgenes_logfoldchange_2.pdf"), height = 3, width = 4)
  
  ## summary of most important genes across clusters
  ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes[DE_Chrimsonposcells_lights_per_cluster_topgenes$gene %in% tail(levels(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene), n=10),], aes(y=gene, x=L1, fill=value))+
  # ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes, aes(y=gene, x=L1, fill=value))+
    scale_fill_gradientn(colours=c(bl,"white", re))+
    geom_tile()+ggtitle('Differentially expressed genes between lights ON/OFF per cell cluster\n(most shared genes)')+
    theme_bw()+labs(x='Cluster', y='Gene')
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes_subsetgenestop10.png"), height = 4, width = 8)
  
  # splitUMAPPlot(Chrimsonposcells, group.by='stim')
  # splitUMAPPlot(Chrimsonposcells, group.by='stim_light')
  
  ## num of DE genes per cluster
  dataset$num_light_DE_genes = table(DE_Chrimsonposcells_lights_per_cluster_topgenes$L1)[match(dataset$seurat_clusters, names(table(DE_Chrimsonposcells_lights_per_cluster_topgenes$L1)))]
  dataset$num_light_DE_genes_any_bool <- !is.na(dataset$num_light_DE_genes)
  dataset$num_light_DE_geneslogp1 <- log(dataset$num_light_DE_genes+1)
  table(dataset$num_light_DE_genes_any_bool)
  
  ## UMAP of lights on/off coloured by feature
  splitUMAPPlot(dataset, group.by='stim_light', colour_by_feature='sr')
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_sr.png"), height = 2.2, width = 4)
  splitUMAPPlot(dataset, group.by='stim_light', colour_by_feature='Hr38')
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_Hr38.png"), height = 2.2, width = 4)
  splitUMAPPlot(dataset, group.by='stim_light', colour_by_feature='Hr38', dim_red = "tsne")
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_Hr38_TSNE.png"), height = 2.2, width = 4)
  splitUMAPPlot(dataset, group.by='stim_light', colour_by_feature='sr', dim_red = "tsne")
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_sr_TSNE.png"), height = 2.2, width = 4)
  
  
  ## percentages of genes which are DE, in each cluster, including all cells, or separating by Chrimson cells
  dataset$num_light_DE_genes[is.na(dataset$num_light_DE_genes)] <- 0
  FeaturePlot(dataset, feature='num_light_DE_genes')+ggtitle('Number of DE genes')
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_numDE_UMAP.png"), height = 2.8, width = 4)
  FeaturePlot(dataset, feature='num_light_DE_geneslogp1')+ggtitle('Number of DE genes')
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_numDElogp1_UMAP.png"), height = 2.8, width = 4)
  FeaturePlot(dataset, feature='num_light_DE_genes_any_bool')+ggtitle('Presence of DE genes')
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_lightONvsOFF_anyDE_UMAP.png"), height = 2.8, width = 4)
  
  # FeaturePlot(Chrimsonposcells, feature='num_light_DE_genes', reduction = 'pca')+ggtitle('Number of DE genes')
  # FeaturePlot(Chrimsonposcells, feature='num_light_DE_genes', reduction = 'tsne')+ggtitle('Number of DE genes')
  
  ## is there a power issue?
  num_DE_genes_cluster = cbind.data.frame(num_DE_genes=as.vector(table(factor(DE_Chrimsonposcells_lights_per_cluster_topgenes$L1, levels=sort(unique(dataset$seurat_clusters))))),
                              cluster_size=as.vector(table(factor(dataset$seurat_clusters, levels=sort(unique(dataset$seurat_clusters))))),
                              cluster=sort(unique(dataset$seurat_clusters)))
  num_DE_genes_cluster
  
  ggplot(num_DE_genes_cluster, aes(x=cluster_size, y=num_DE_genes))+
    geom_point()+geom_smooth(method = "lm")+theme_bw()+labs(x='Size of cluster', y='Number of DE genes in cluster')
  ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/", dataset_name, "/", dataset_name, "_cor_numDEgenes_clustersize.png"), height = 3, width = 3)

}

give_DE_analysis_lightON_lightOFF(Chrimsonposcells, 'chrimsonpos')
give_DE_analysis_lightON_lightOFF(dataset = Chrimsonnegcells, dataset_name = 'chrimsonneg')
give_DE_analysis_lightON_lightOFF(dataset = combined_dataset, dataset_name = 'allcells')

add_genenames <- function(i){
  lapply(i, function(j) cbind.data.frame(gene=rownames(j), j))
}
DE_lights_per_cluster <- lapply(c('chrimsonpos', 'chrimsonneg', 'allcells'), function(i){
  x <- readRDS(paste0("~/projects/lighting/data/robjects/", input_objs, "/", i, "_DEgenes_per_cluster.RDS"))
  add_genenames(x)
})
names(DE_lights_per_cluster) <- c('chrimsonpos', 'chrimsonneg', 'allcells')

# DE_lights_per_cluster
head(melt(DE_lights_per_cluster$chrimsonneg, id.vars=c('gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj')))
DE_lights_per_cluster_molten <- (melt(DE_lights_per_cluster, id.vars=c('gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj')))

## summaries
sort(table(DE_lights_per_cluster_molten[(DE_lights_per_cluster_molten$L1 == "chrimsonpos") & (DE_lights_per_cluster_molten$p_val_adj <= 0.05),'L2']))
median(table(DE_lights_per_cluster_molten[(DE_lights_per_cluster_molten$L1 == "chrimsonpos") & (DE_lights_per_cluster_molten$p_val_adj <= 0.05),'L2']))
sort(table(DE_lights_per_cluster_molten[(DE_lights_per_cluster_molten$L1 == "chrimsonpos") & (DE_lights_per_cluster_molten$p_val_adj <= 0.05),'gene']))

DE_lights_per_cluster_molten_top_common <- DE_lights_per_cluster_molten[DE_lights_per_cluster_molten$gene %in% names(tail(sort(table(DE_lights_per_cluster_molten$gene)), n=20)),]

DE_lights_per_cluster_molten_top_common$L1L2 <- apply(DE_lights_per_cluster_molten_top_common[,c('L1', 'L2')],1, paste0, collapse='-')
DE_lights_per_cluster_molten_top_common$signif <- DE_lights_per_cluster_molten_top_common$p_val_adj <= 0.05
DE_lights_per_cluster_molten_top_common_dcast <- (dcast(DE_lights_per_cluster_molten_top_common[,c('L1L2', 'gene', 'signif')], gene~L1L2, value.var = "signif"))
rownames(DE_lights_per_cluster_molten_top_common_dcast) <- DE_lights_per_cluster_molten_top_common_dcast$gene; DE_lights_per_cluster_molten_top_common_dcast[,1] <- NULL
DE_lights_per_cluster_molten_top_common_dcast[is.na(DE_lights_per_cluster_molten_top_common_dcast)] <- FALSE

DE_lights_per_cluster_molten_top_common_pheatmap <- pheatmap::pheatmap(apply(DE_lights_per_cluster_molten_top_common_dcast, 1, as.numeric))
(DE_lights_per_cluster_molten_top_common_pheatmap)

DE_lights_per_cluster_molten_top_common_labels <- colnames(DE_lights_per_cluster_molten_top_common_dcast)[DE_lights_per_cluster_molten_top_common_pheatmap$tree_row$order]
DE_lights_per_cluster_molten_top_common_labelsgenes <- rownames(DE_lights_per_cluster_molten_top_common_dcast)[DE_lights_per_cluster_molten_top_common_pheatmap$tree_col$order]

relevel_groups <- function(i){
  if(i == 'chrimsonneg'){
    'Chrimson negative'
  }else if(i == 'chrimsonpos'){
    'Chrimson positive'
  }else if(i == 'allcells'){
    'All cells'
  }
}
DE_lights_per_cluster_molten_top_common$L1 <- sapply(DE_lights_per_cluster_molten_top_common$L1, relevel_groups)
ggplot(DE_lights_per_cluster_molten_top_common, aes(x=factor(gene, levels=DE_lights_per_cluster_molten_top_common_labelsgenes),
                                                    y=factor(L1L2, levels=DE_lights_per_cluster_molten_top_common_labels),
                                                    fill=avg_log2FC))+
  geom_tile()+facet_wrap(.~L1, ncol=1, scales = "free_y")+theme_bw()+
  jcolors::scale_fill_jcolors_contin(palette = "pal3")+
  labs(x='Gene', y='Cluster')+
  theme(axis.text.x = element_text(angle = 30, hjust=1))+
  theme(#axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(legend.position = "bottom")
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/lightONvsOFF_summarygenes_tile.pdf"), height = 5, width = 7)

ggplot(DE_lights_per_cluster_molten_top_common[DE_lights_per_cluster_molten_top_common$L1 %in% c('Chrimson negative', 'Chrimson positive'),],
       aes(y=factor(gene, levels=DE_lights_per_cluster_molten_top_common %>% group_by(gene) %>% dplyr::summarise(median(avg_log2FC)) %>% dplyr::arrange(`median(avg_log2FC)`) %>% dplyr::select(gene) %>% unlist),
                                                    x=avg_log2FC, fill=L1, alpha=0.2), lty=3)+geom_joy()+theme_bw()+geom_vline(xintercept = 0)+
  scale_fill_jcolors(palette = "pal6")+labs(x='Average log-fold change across cluster', y='Gene', fill='')+
  theme(legend.position = "bottom")+guides(alpha="none")
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/lightONvsOFF_summarygenes_joy.pdf"), height = 5, width = 7)

dcastout_to_df <- function(i){
  rownames(i) <- i[,1]
  i[,-1]
}

add_average_median_SD <- function(i){
  a <- cbind(median=dcast(i, gene~L1, value.var = 'med'), mean=dcast(i, gene~L1, value.var = 'mean')[,-1],
        sd=dcast(i, gene~L1, value.var = 'sd')[,-1])
  a$col = as.factor(as.numeric(a[,'mean.Chrimson positive',] > a[,'mean.Chrimson negative']))
  a
}
head(add_average_median_SD(DE_lights_per_cluster_molten_top_common %>% filter( (L1 != 'All cells') & signif ) %>% group_by(gene, L1) %>%
                        dplyr::summarise(mean=mean(avg_log2FC), med=median(avg_log2FC), sd=sd(avg_log2FC))))

give_errorbar_plot_LFC <- function(i){
  ggplot(i)+
  geom_point(aes(x=`median.Chrimson negative`, y=`median.Chrimson positive`), col='blue')+
  geom_point(aes(x=`mean.Chrimson negative`, y=`mean.Chrimson positive`))+
  geom_errorbar(aes(x=`mean.Chrimson negative`,
                    ymin=`mean.Chrimson positive`-`sd.Chrimson positive`,
                    ymax=`mean.Chrimson positive`+`sd.Chrimson positive`), alpha=0.5)+
  geom_errorbar(aes(y=`mean.Chrimson positive`,
                    xmin=`mean.Chrimson negative`-`sd.Chrimson negative`,
                    xmax=`mean.Chrimson negative`+`sd.Chrimson negative`), alpha=0.5)+
  geom_segment(aes(x=`mean.Chrimson negative`, y=`mean.Chrimson positive`, xend=`median.Chrimson negative`, yend=`median.Chrimson positive`))+
  geom_abline(intercept = 0, slope = 1, lty='dashed')+theme_bw()+
  geom_label_repel(aes(x=`mean.Chrimson negative`, y=`mean.Chrimson positive`, label=median.gene,
                       col=col
                       ))+
  scale_color_manual(values=c('1'='black', '0'='red'))+
    guides(col='none')+
    labs(x='LFC of GE in Chrimson- clusters',
         y='LFC of GE in Chrimson+ clusters')
}

give_errorbar_plot_LFC(add_average_median_SD(DE_lights_per_cluster_molten_top_common %>% filter( (L1 != 'All cells') ) %>% group_by(gene, L1) %>%
                                               dplyr::summarise(mean=mean(avg_log2FC), med=median(avg_log2FC), sd=sd(avg_log2FC))))
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/lightONvsOFF_chrimsonpos_chrimsonneg_LFCcomparison.pdf"), height = 3, width = 3)

give_errorbar_plot_LFC(add_average_median_SD(DE_lights_per_cluster_molten_top_common %>% filter( (L1 != 'All cells') & signif ) %>% group_by(gene, L1) %>%
                         dplyr::summarise(mean=mean(avg_log2FC), med=median(avg_log2FC), sd=sd(avg_log2FC))))
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/lightONvsOFF_chrimsonpos_chrimsonneg_LFCcomparison_onlyDEclusters.pdf"), height = 3, width = 3)

DE_lights_per_cluster_molten_top_commonnumDE <- (DE_lights_per_cluster_molten_top_common) %>% dplyr::filter(p_val_adj <= 0.05) %>% dplyr::group_by(L1L2, L1) %>% dplyr::count((L1L2))
DE_lights_per_cluster_molten_top_commonnumclusters <- (DE_lights_per_cluster_molten_top_common) %>% dplyr::filter(p_val_adj <= 0.05) %>% dplyr::group_by(L1) %>% dplyr::count((gene))
DE_lights_per_cluster_molten_top_commonnumDE

colnames(DE_lights_per_cluster_molten_top_commonnumDE) <- c('L1L2', 'L1', 'cluster', 'n')
ggplot(DE_lights_per_cluster_molten_top_commonnumDE, aes(x=factor(cluster, levels=DE_lights_per_cluster_molten_top_commonnumDE$cluster[order(DE_lights_per_cluster_molten_top_commonnumDE$n)]), y=n))+geom_bar(stat = "identity")+
  facet_wrap(.~L1, ncol=1)+
  theme(axis.text.x=element_blank(),
    axis.ticks.x=element_blank())+labs(x='Clusters', y='Number of DE genes')
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/DE_per_cluster.pdf"), height = 5, width = 4)

DE_lights_per_cluster_molten_top_commonnumclusterssum <- DE_lights_per_cluster_molten_top_commonnumclusters %>% group_by(`(gene)`) %>% summarise(sum(n))
num_cells_datasets <- data.frame(dataset=c( "All cells", "Chrimson negative", "Chrimson positive"),
                                 numcells=c(dim(combined_dataset@assays$RNA)[2], dim(Chrimsonnegcells@assays$RNA)[2], dim(Chrimsonposcells@assays$RNA)[2]),
                                 numclusters=c(length(levels(combined_dataset$seurat_clusters)), length(levels(Chrimsonnegcells$seurat_clusters)), length(levels(Chrimsonposcells$seurat_clusters))))
DE_lights_per_cluster_molten_top_commonnumclusters$totalcells <- num_cells_datasets$numcells[match(DE_lights_per_cluster_molten_top_commonnumclusters$L1, num_cells_datasets$dataset)]
DE_lights_per_cluster_molten_top_commonnumclusters$totalclusters <- num_cells_datasets$numclusters[match(DE_lights_per_cluster_molten_top_commonnumclusters$L1, num_cells_datasets$dataset)]
DE_lights_per_cluster_molten_top_commonnumclusters$percent=DE_lights_per_cluster_molten_top_commonnumclusters$n/DE_lights_per_cluster_molten_top_commonnumclusters$totalclusters*100

ggplot(DE_lights_per_cluster_molten_top_commonnumclusters, aes(x=factor(`(gene)`, levels=DE_lights_per_cluster_molten_top_commonnumclusterssum$`(gene)`[order(DE_lights_per_cluster_molten_top_commonnumclusterssum$`sum(n)`)]), y=n))+geom_bar(stat = "identity")+
  facet_wrap(.~L1, ncol=1)+
  labs(x='Genes', y='Number of DE clusters where the gene is DE')+theme_bw()+
  theme(axis.text.x = element_text(angle = 30, hjust=1))
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/DE_per_gene.pdf"), height = 5, width = 4)

ggplot(DE_lights_per_cluster_molten_top_commonnumclusters, aes(x=factor(`(gene)`, levels=DE_lights_per_cluster_molten_top_commonnumclusterssum$`(gene)`[order(DE_lights_per_cluster_molten_top_commonnumclusterssum$`sum(n)`)]),
                                                               y=percent))+geom_bar(stat = "identity")+
  facet_wrap(.~L1, ncol=1)+
  labs(x='Genes', y='Percentage of DE clusters where the gene is DE (%)')+theme_bw()+
  theme(axis.text.x = element_text(angle = 30, hjust=1))
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/DE_per_gene_percent.pdf"), height = 5, width = 4)


## check, for each cluster in the all-cells dataset, to what extent it is composed of chrimson+cells
chrimson_exprs_per_cluster0 <- sapply(levels(combined_dataset$seurat_clusters), function(cluster_it){
  mean(combined_dataset@assays$RNA['Chrimson'][,combined_dataset$seurat_clusters == cluster_it] > 0)
})

hist(chrimson_exprs_per_cluster0)

chrimson_exprs_per_cluster <- cbind.data.frame(chrimson_exprs_per_cluster0, cluster=names(chrimson_exprs_per_cluster0), ordered=as.integer(rank(-chrimson_exprs_per_cluster0)))
chrimson_exprs_per_cluster <- chrimson_exprs_per_cluster[order(chrimson_exprs_per_cluster$chrimson_exprs_per_cluster0),]
plot(chrimson_exprs_per_cluster$ordered, chrimson_exprs_per_cluster$chrimson_exprs_per_cluster0)
ggplot(chrimson_exprs_per_cluster,
       aes(x = ordered, label=cluster,
           # x=factor(cluster, levels=cluster[order(chrimson_exprs_per_cluster, decreasing =T)]),
           y=chrimson_exprs_per_cluster0, group=1))+
  geom_point()+geom_line()+theme_bw()+lims(x=c(0,10))+geom_label_repel()+
  labs(x='Cluster', y='Fraction of cells in cluster expressing Chrimson')
ggsave(paste0("~/projects/lighting/3_results/", input_objs, "/fraction_chrimson_in_clusters.pdf"), height = 4, width = 4)

combined_dataset$fraction_Chrimson_in_cluster

# save.image(file = "~/projects/lighting/data/robjects/image_DEA_Charly_DE.RData")
# load(file = "~/projects/lighting/data/robjects/image_DEA_Charly_DE.RData")
# .rs.restartR()

## --------------------------------------------------------- ##
## Analyses without the transposable elements
# load(file = "~/projects/lighting/data/robjects/image_DEA_Charly.RData")

if(input_objs == 'CharlyCellRanger'){
  
  stop("Paths to save files/images have not been changed. Change them before running this code\n")
  
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
  
}

system(paste0("ls ~/projects/lighting/data/robjects/", input_objs, "/image.RData"))
save.image(file = paste0("~/projects/lighting/data/robjects/", input_objs, "/image.RData"))

