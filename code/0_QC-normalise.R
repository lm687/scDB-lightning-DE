##https://pachterlab.github.io/kallistobustools/tutorials/kb_getting_started/R/kb_intro_2_R/

rm(list = ls())
set.seed(1234)
local <- F

##--------------------------------------------------------------------------------------------------##
library(optparse)
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
if(local){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
  
  opt <- list()
  # opt$genome <- 'KB-dmel649ChrimsonV2'
  # opt$folder_alignment_files <- '/home/l/lmorrill/projects/lighting/data/1_alignment/KallistoOUT/dmel649ChrimsonV2/'
  
  opt$genome <- 'KB-dmel649Chrimson'
  opt$folder_alignment_files <- '/home/l/lmorrill/projects/lighting/data/1_alignment/KallistoOUT/dmel649Chrimson/'
  
  # opt$folder_alignment_files <- '/home/l/lmorrill/projects/lighting/data/1_alignment/cellrangerOUT/dmel649Chrimson/' ## does not have Chrimson
  # opt$genome <- 'CR-dmel649Chrimson'
  # opt$folder_alignment_files <- '/home/l/lmorrill/projects/lighting/data/1_alignment/cellrangerOUT/dmel649CORRECTED/' ## does not have Chrimson
}else{
  option_list = list(
    make_option(c("--genome"), type="character", default="None",
                help="name of alignment files (corresponding to the genome used for alignment)", metavar="character"),
    make_option(c("--folder_alignment_files"), type="character", default="None",
                help="name of folder which contains all the alignment files, each in a subfolder", metavar="character")
  )
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
}
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
folder_results <- paste0("~/projects/lighting/results/", opt$genome, '/') ## CCB cluster
robject_folder <- paste0("~/projects/lighting/robjects/", opt$genome, '/') ## CCB cluster

system(paste0("mkdir -p ", folder_results))
system(paste0("mkdir -p ", robject_folder))

folder_results <- paste0(folder_results, "0_QC-normalise", '/') ## CCB cluster
robject_folder <- paste0(robject_folder, "0_QC-normalise", '/') ## CCB cluster
system(paste0("mkdir -p ", folder_results))
system(paste0("mkdir -p ", robject_folder))

##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
library(DropletUtils)
library(Matrix)
library(tidyverse)
library(Seurat)
library(ggpointdensity)
library(scico)
library(scales)
library(gridExtra)
library(jcolors)
library(ggrepel)
library(fgsea)
library(reshape2)
library(celda) ## decontX
library(SoupX)
library(GGally)
library(scDblFinder)
theme_set(theme_bw())
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
threshold_ambient_soupx = .3
threshold_ambient_decontx = .3
threshold_upper_nCount = 10000
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## functions
source("helper_functions.R")
source("helper_kallisto.R")
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
type_alignment <- if(grepl('^Kallisto', basename(dirname(opt$folder_alignment_files)))){
  'Kallisto'
}else if(grepl('^cellranger|Cellranger', basename(dirname(opt$folder_alignment_files)))){
  'CellRanger'
}else{
  stop('Unrecognised alignment folder\n')
}
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
samples <- apply(expand.grid(paste0('G', 1:4), paste0('_rep', 1:2)), 1, paste0, collapse='')
samples <- gtools:::mixedsort(samples)
samples
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## read counts; check presence of chrimson
if(type_alignment == 'Kallisto'){
  res_mat <- lapply(samples, function(i){
    read_count_output(paste0(opt$folder_alignment_files, "/",
                             i, "/counts_unfiltered"), name = "cells_x_genes")
  })
  names(res_mat) <- samples
  res_mat0 <- res_mat ## no filter. Used for SoupX
}else if(type_alignment == 'CellRanger'){
  res_mat <- lapply(samples, function(i){
    Read10X(paste0(opt$folder_alignment_files, i, '/', i, "/outs/filtered_feature_bc_matrix"))
  })
  res_mat0 <- lapply(samples, function(i){
    Read10X(paste0(opt$folder_alignment_files, i, '/', i, "/outs/raw_feature_bc_matrix"))
  })
  
}else{
  stop('Unrecognised <type_alignment>\n')
}
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
chrimson_features_cells <- sapply(res_mat, function(i) table(i['FBto0000555',]))
stopifnot(all(sapply(res_mat, function(i) sum(i['FBto0000555',]) > 0))) ## stop if there aren't chrimson counts in all samples
chrimson_features_cells_percent <- sapply(chrimson_features_cells, function(i) 100*(1-(i[1]/sum(i)))) ## percentage of chrimson pos. cells
## these results are the same as those with the dmel649Chrimson reference genome
chrimson_features_cells_percent*100
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## filtering empty droplets

if(type_alignment == 'Kallisto'){
  tot_counts <- lapply(res_mat, Matrix::colSums)
  bc_rank <- lapply(res_mat, barcodeRanks, lower = 10)
  names(tot_counts) <- samples
  names(bc_rank) <- samples
  
  ## knee plots
  options(repr.plot.width=9, repr.plot.height=6)
  knee_plots <- lapply(samples, function(i) knee_plot(bc_rank[[i]])+ggtitle(i))
  names(knee_plots) <- samples
  knee_plots <- knee_plots[gtools::mixedorder(names(knee_plots))]
  # do.call('grid.arrange', knee_plots)
  
  inflections <- sapply(samples, function(i) metadata(bc_rank[[i]])$inflection)
  for(i in 1:length(knee_plots)){
    knee_plots[[i]] <- knee_plots[[i]]+geom_vline(xintercept = inflections[i], col='red')
  }
  
  pdf(paste0(folder_results, "/1_kneeplots.pdf"), height = 5.5, width = 3.5)
  grid.arrange(grobs=knee_plots, nrow=4)
  dev.off()
  
  ## filtering
  res_mat <- lapply(samples, function(i) res_mat[[i]][, tot_counts[[i]] > metadata(bc_rank[[i]])$inflection])
  res_mat <- lapply(res_mat, function(i) i[Matrix::rowSums(i) > 0,])
  names(res_mat) <- samples
  sapply(res_mat, dim)
}

##--------------------------------------------------------------------------------------------------##
## Rename genes, if necessary
name_conversion_file <- give_name_conversion_file()
name_conversion_file <- rbind(name_conversion_file, c('FBto0000555', 'Chrimson')) ## Adding Chrimson
names(name_conversion_file) <- c('gene', 'gene_symbol')
head(name_conversion_file)

for(i in samples){
  rownames(res_mat[[i]]) <- name_conversion_file$gene_symbol[match(rownames(res_mat[[i]]), name_conversion_file$gene)]
}

for(i in samples){
  rownames(res_mat0[[i]]) <- name_conversion_file$gene_symbol[match(rownames(res_mat0[[i]]), name_conversion_file$gene)]
  stopifnot(!any(is.na(rownames(res_mat0[[i]]))))
}

##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Create seurat object (for now, one per sample)
seu <- lapply(res_mat, function(i) CreateSeuratObject(i, min.cells = 3, min.features = 200))
for(i in 1:length(seu)){ ## add name
  seu[[i]]$orig.ident <- samples[i]
}
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Add percent mt
for(i in samples){
  seu[[i]][["percent.mt"]] <- PercentageFeatureSet(seu[[i]], pattern = "^mt:")
}

sapply(samples, function(i) mean(seu[[i]]$percent.mt, na.rm=T))
ggplot(rownames_as_col(melt((sapply(samples, function(i) mean(seu[[i]][["percent.mt"]][,1], na.rm=T))))),
       aes(x=name, y=value, group=1))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle = 30, hjust=1))+labs(x='Sample', y='Percentage of mt')
ggsave(paste0(folder_results, "/1_percent_mt.pdf"), height = 3.5, width = 3.5)

# vlplots <- lapply(seu, function(i) VlnPlot(i, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1))

ggplot(melt(lapply(seu, function(i) i[[c("nFeature_RNA", "nCount_RNA", "percent.mt")]])), aes(x=L1, y=value))+
  geom_violin()+facet_wrap(.~variable, scales='free', nrow=3)+scale_y_log10()+
  theme(axis.text.x = element_text(angle = 30, hjust=1))+labs(x='Sample')
ggsave(paste0(folder_results, "/1_violinplots_QC.pdf"), height = 5.5, width = 3.5)

options(repr.plot.width=9, repr.plot.height=6)
plots_UMI_genesdetected <- lapply(samples, function(i) ggplot(seu[[i]]@meta.data, aes(nCount_RNA, nFeature_RNA)) +
                                    geom_hex(bins = 100) +
                                    scale_fill_scico(palette = "devon", direction = -1, end = 0.9) +
                                    scale_x_log10(breaks = breaks_log(12)) + 
                                    scale_y_log10(breaks = breaks_log(12)) + annotation_logticks() +
                                    labs(x = "Total UMI counts", y = "Number of genes detected") +
                                    theme(panel.grid.minor = element_blank())+ggtitle(i))
names(plots_UMI_genesdetected) <- samples

plots_UMI_genesdetected <- plots_UMI_genesdetected[gtools::mixedorder(names(plots_UMI_genesdetected))]
pdf(paste0(folder_results, "/1_UMI_ngenes.pdf"), height = 7.5, width = 7.5)
grid.arrange(grobs=plots_UMI_genesdetected, ncol=2)
dev.off()


ggplot(do.call('rbind', lapply(seu, function(i) i@meta.data)), aes(nCount_RNA, nFeature_RNA, col=orig.ident)) +
  geom_point() +
  scale_x_log10(breaks = breaks_log(12)) + 
  scale_y_log10(breaks = breaks_log(12)) + annotation_logticks() +
  labs(x = "Total UMI counts", y = "Number of genes detected") +
  theme(panel.grid.minor = element_blank())+ggtitle('All samples')+geom_smooth()
ggsave(paste0(folder_results, "/1_UMI_ngenes_2.pdf"), height = 4.5, width = 4.5)

plots_UMI_mt <- lapply(samples, function(i) ggplot(seu[[i]]@meta.data, aes(nCount_RNA, percent.mt)) +
                         geom_pointdensity() +
                         scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
                         labs(x = "Total UMI counts", y = "Percentage mitochondrial")+ggtitle(i))

pdf(paste0(folder_results, "/1_UMI_mt.pdf"), height = 7.5, width = 7.5)
grid.arrange(grobs=plots_UMI_mt, ncol=2)
dev.off()

plots_UMI_mt_v2 <- lapply(samples, function(i) ggplot(seu[[i]]@meta.data, aes(nCount_RNA, percent.mt)) +
                            geom_pointdensity() +
                            scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
                            #scale_x_log10()+scale_y_log10()+
                            lims(x=c(min(seu[[i]]@meta.data$nCount_RNA),500), y=c(0,10))+
                            labs(x = "Total UMI counts", y = "Percentage mitochondrial")+ggtitle(i))

pdf(paste0(folder_results, "/1_UMI_mt_v2.pdf"), height = 7.5, width = 7.5)
grid.arrange(grobs=plots_UMI_mt_v2, ncol=2)
dev.off()

# We filter cells with more than 3% mitochondrial content based on the plot above.
# I use a shared value of 3% mt DNA
seu <- lapply(seu, subset, subset = percent.mt < 3)

##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Centering and scaling data matrix
seu <- lapply(seu, function(i) NormalizeData(i) %>% ScaleData())
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
# Identify highly variable genes

## With vst (by default)
seu <- lapply(seu, FindVariableFeatures, nfeatures = 2000)
top10 <-  lapply(seu, function(i) head(VariableFeatures(i), n=10))
plot1 <-  lapply(seu, VariableFeaturePlot, log = FALSE)
HVG_plot <- lapply(samples, function(i) LabelPoints(plot1[[i]], points = top10[[i]], repel = TRUE))

# do.call('grid.arrange', HVG_plot)

ggplot(melt(!apply(firstcol_to_rownames(dcast(melt(top10), value~L1)), 2, is.na)) %>% dplyr::filter(value),
       aes(y=Var1, x=Var2, fill=value))+geom_tile()+labs(fill='Is Top10 HVG')+
  theme(axis.text.x = element_text(angle = 30, hjust=1))+
  scale_fill_manual(values = c("darksalmon"))+labs(x='Sample', y='HVG')+theme(legend.position = "bottom")
ggsave(paste0(folder_results, "/1_HVG_top10_tile.pdf"), height = 6.5, width = 4.5)

pdf(paste0(folder_results, "/1_HVG_top10_plot.pdf"), height = 10.5, width = 10.5)
grid.arrange(grobs=lapply(samples, function(i) VariableFeaturePlot(seu[[i]])+ggtitle(i)), ncol=2)
dev.off()

## With MVP
seu <- lapply(seu, FindVariableFeatures, nfeatures = 2000, selection.method='mvp')
top10_MVP <-  lapply(seu, function(i) head(VariableFeatures(i), n=10))
ggplot(melt(!apply(firstcol_to_rownames(dcast(melt(top10_MVP), value~L1)), 2, is.na)) %>% dplyr::filter(value),
       aes(y=Var1, x=Var2, fill=value))+geom_tile()+labs(fill='Is Top10 HVG (MVP)')+
  theme(axis.text.x = element_text(angle = 30, hjust=1))+
  scale_fill_manual(values = c("salmon"))+labs(x='Sample', y='HVG (MVP)')+theme(legend.position = "bottom")
ggsave(paste0(folder_results, "/1_HVG_MVP_top10_tile.pdf"), height = 6.5, width = 4.5)

pdf(paste0(folder_results, "/1_HVG_MVP_top10_plot.pdf"), height = 10.5, width = 10.5)
grid.arrange(grobs=lapply(samples, function(i) VariableFeaturePlot(seu[[i]])+ggtitle(i)), ncol=2)
dev.off()

## Using MVP from now on

##--------------------------------------------------------------------------------------------------##
## PCA

plot_elbow_vlines <- function(seurat_obj){
  pca_obj <- Reductions(seu$G1_rep1, 'pca')
  ggplot(cbind.data.frame(x=1:length(pca_obj@stdev), y=pca_obj@stdev))+
    geom_segment(aes(x=x,xend=x,y=0,yend=y))+geom_point(aes(x=x, y=y))+labs(x='PC', y='Standard deviations')
}

seu <- lapply(seu, function(i) RunPCA(i, verbose = FALSE, npcs = 50)) # uses HVG by default

elbows <- lapply(samples, function(i) plot_elbow_vlines(seu[[i]])+ggtitle(i))

pdf(paste0(folder_results, "/1_pca_elbow.pdf"), height = 4.5, width = 6.5)
grid.arrange(grobs=elbows, ncol=4)
dev.off()

pcas <- lapply(samples, function(i) PCAPlot(seu[[i]])+ggtitle(i))

pdf(paste0(folder_results, input_objs, "/1_pcas.pdf"), height = 15, width = 7.5)
grid.arrange(grobs=pcas, ncol=2)
dev.off()
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Ambient and doublets
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## DecontX
sce <- lapply(seu, function(i) SingleCellExperiment(assays=list(counts=i@assays$RNA@counts)))
sce <- lapply(sce, scuttle::logNormCounts)

## Ambient
decontx_res <- lapply(sce, celda::decontX)
decontx_res
umap_decontX <- lapply(decontx_res, reducedDim, "decontX_UMAP")

pdf(paste0(folder_results, "/2_decontX_UMAP.pdf"), height = 7.5, width = 4.5)
grid.arrange(grobs= lapply(samples, function(i){
  plotDimReduceCluster(x = decontx_res[[i]]$decontX_clusters, 
                       dim1 = umap_decontX[[i]][, 1], dim2 = umap_decontX[[i]][, 2])+ggtitle(i)
}), ncol=2)
dev.off()

pdf(paste0(folder_results, "/2_decontX_UMAP_contamination.pdf"), height = 10.5, width = 7.5)
grid.arrange(grobs= lapply(samples, function(i){
  plotDecontXContamination(x = decontx_res[[i]])+ggtitle(i)
}), ncol=2)
dev.off()
## really high contamination

##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## DecontX with clusters
## find neighbours/clustering
seu <- lapply(seu, FindNeighbors, dims = 1:30)
seu <- lapply(seu, FindClusters)
seu <- lapply(seu, RunUMAP, dims = 1:30)

decontx_res_clusters <- lapply(samples, function(i) celda::decontX(x = sce[[i]],
                                                                   z=seu[[i]]$seurat_clusters))
umap_decontX_clusters <- lapply(decontx_res_clusters, reducedDim, "decontX_UMAP")
names(decontx_res_clusters) <- names(umap_decontX_clusters) <- samples

pdf(paste0(folder_results, "/2_decontXclusters_UMAP.pdf"), height = 7.5, width = 4.5)
grid.arrange(grobs= lapply(samples, function(i){
  plotDimReduceCluster(x = decontx_res_clusters[[i]]$decontX_clusters, 
                       dim1 = umap_decontX_clusters[[i]][, 1], 
                       dim2 = umap_decontX_clusters[[i]][, 2])+ggtitle(i)
}), ncol=2)
dev.off()

pdf(paste0(folder_results, "/2_decontXclusters_UMAP_contamination.pdf"), height = 10.5, width = 7.5)
grid.arrange(grobs= lapply(samples, function(i){
  plotDecontXContamination(x = decontx_res[[i]])+ggtitle(i)
}), ncol=2)
dev.off()

##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## soupX
# soupx <- SoupX::autoEstCont(sce$G1_rep1)
sc = sapply(samples, function(i) SoupX::SoupChannel(res_mat0[[i]][rownames(seu[[i]]@assays$RNA@counts),],
                                                    seu[[i]]@assays$RNA@counts),
            simplify = F, USE.NAMES = T)
sc = sapply(samples, function(i) setClusters(sc[[i]], clusters = seu[[i]]$seurat_clusters),
            simplify = F, USE.NAMES = T)
sc = sapply(sc, autoEstCont, simplify = F, USE.NAMES = T)
soupx_out = sapply(sc, adjustCounts, simplify = F, USE.NAMES = T)

##--------------------------------------------------------------------------------------------------##


##--------------------------------------------------------------------------------------------------##
## add this information to the seurat objects

## add DecontX and SoupX ambient estimates in the metadata
for(i in samples){
  seu[[i]]$decontX_contamination = decontx_res[[i]]$decontX_contamination
  seu[[i]]$decontX_clusters = decontx_res[[i]]$decontX_clusters
}

for(i in samples){
  seu[[i]]$decontXclusters_contamination = decontx_res_clusters[[i]]$decontX_contamination
  seu[[i]]$decontXclusters_clusters = decontx_res_clusters[[i]]$decontX_clusters
}

## summarise soupX ambient: get fraction of ambient
for(i in samples){
  stopifnot(all(colnames(seu[[i]]@assays$RNA@counts) == colnames(soupx_out[[i]])))
  seu[[i]]$soupX_frac_ambient = (1- colSums(soupx_out[[i]])/colSums(seu[[i]]@assays$RNA@counts))
}


## add SoupX as a new assay
for(i in samples){
  seu[[i]][["SoupX_adjusted"]] <- CreateAssayObject(data = soupx_out[[i]])
}

for(i in samples){
  cat(Assays(seu[[i]]), '\n')
}


sc = sapply(samples, function(i) setDR(sc[[i]], DR=seu[[i]]@reductions$umap@cell.embeddings),
            simplify = F, USE.NAMES = T)

pdf(paste0(folder_results, "/2_soupX_Vmat_MarkerMap.pdf"), height = 8, width = 8)
for(i in samples){
  print(plotMarkerMap(sc[[i]], "Vmat"))
}
dev.off()
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Comparison of SoupX and DecontX estimates of ambient RNA contribution
pdf(paste0(folder_results, "/2_soupX_DecontX_UMAP_contamination.pdf"), height = 3.5, width = 10)
for(i in samples){
  print(FeaturePlot(object = seu[[i]], features = c('decontX_contamination',
                                                    'decontXclusters_contamination',
                                                    'soupX_frac_ambient'), ncol=3)+labs(x=paste0('UMAP_1 ', i)))
}
dev.off()

pdf(paste0(folder_results, "/2_soupX_DecontX_ggpairs_contamination.pdf"), height = 8, width = 8)
for(i in samples){
  print(ggpairs(seu[[i]]@meta.data[,c('decontX_contamination',
                                      'decontXclusters_contamination',
                                      'soupX_frac_ambient')]))
}
dev.off()
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Doublets
sce_scdblfinder <- sapply(samples, function(i){
  scDblFinder(sce = sce[[i]])
}, simplify = F, USE.NAMES = T)

## adding this information to the seurat object
for(i in samples){
  stopifnot(all(colnames(seu[[i]]@assays$RNA@counts) == colnames(sce_scdblfinder[[i]])))
  seu[[i]]$scDblFinder_class = sce_scdblfinder[[i]]$scDblFinder.class
}

normalise_cl(sapply(samples, function(i) table(seu[[i]]$scDblFinder_class)))

## plot doublets

pdf(paste0(folder_results, "/2_doublets_UMAP.pdf"), height = 3.5, width = 10)
for(i in samples){
  Idents(seu[[i]]) <- 'scDblFinder_class'
  print(UMAPPlot(object = seu[[i]])+ggtitle(i))
}
dev.off()

##--------------------------------------------------------------------------------------------------##


##--------------------------------------------------------------------------------------------------##
## Filter cells based on
#' 1. Doublet classification
#' 2. Ambient (both SoupX and DecontX)
#' 3. Exclude droplets with extreme nCount values

for(i in samples){
  #' 1. Doublet classification
  seu[[i]] <- seu[[i]][,seu[[i]]$scDblFinder_class == 'singlet']
  #' 2. Ambient (both SoupX and DecontX)
  seu[[i]] <- seu[[i]][,seu[[i]]$soupX_frac_ambient < threshold_ambient_soupx]
  seu[[i]] <- seu[[i]][,seu[[i]]$decontXclusters_contamination < threshold_ambient_decontx]
  #' 3. Exclude droplets with extreme nCount values
  seu[[i]] <- seu[[i]][,seu[[i]]$nCount_RNA <= threshold_upper_nCount]
}
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Re-normalise data

## Using Lognormalisation
## Recompute PCA

for(i in samples){
  DefaultAssay(seu[[i]]) <- 'RNA'
}
seu <- lapply(seu, FindVariableFeatures, nfeatures = 2000, selection.method='mvp')
seu <- lapply(seu, function(i) NormalizeData(i) %>% ScaleData())
seu <- lapply(seu, function(i) RunPCA(i, verbose = FALSE, npcs = 30, assay = 'RNA',
                                      reduction.key = 'PCA', reduction.name='pca')) # uses HVG by default

## Using SCT
seu <- lapply(seu, SCTransform)
sapply(seu, DefaultAssay)
seu <- lapply(seu, function(i) RunPCA(i, verbose = FALSE, npcs = 30, assay = 'SCT',
                                      reduction.key = 'SCTPCA', reduction.name='SCTpca')) # uses HVG by default

names(seu[[1]]@reductions) ## UMAP is the previous

## Recompute UMAP
seu <- lapply(seu, function(i) RunUMAP(i, verbose = FALSE, dims = 1:30,
                                       reduction = 'pca', reduction.key = 'UMAP_',
                                       reduction.name='umap'))
seu <- lapply(seu, function(i) RunUMAP(i, verbose = FALSE, dims = 1:30,
                                       reduction = 'SCTpca', reduction.key = 'SCTUMAP_',
                                       reduction.name='SCTumap'))

names(seu[[1]]@reductions) ## four reductions

##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Re-plot UMAPs
pdf(paste0(folder_results, "/3_filtered_UMAP_nCount.pdf"), height = 3.5, width = 10)
for(i in samples){
  print(FeaturePlot(object = seu[[i]], features = 'nCount_RNA', reduction = 'umap')+ggtitle(i))
}
dev.off()

pdf(paste0(folder_results, "/3_filtered_sctUMAP_nCount.pdf"), height = 3.5, width = 10)
for(i in samples){
  print(FeaturePlot(object = seu[[i]], features = 'nCount_RNA', reduction = 'SCTumap')+ggtitle(i))
}
dev.off()
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
for(i in samples){
  DefaultAssay(seu[[i]]) <- 'RNA'
}
seu <- lapply(seu, function(i) i %>% FindNeighbors %>% FindClusters )

pdf(paste0(folder_results, "/3_filtered_violin_nCount.pdf"), height = 3.5, width = 7)
for(i in samples){
  print(VlnPlot(object = seu[[i]], features = 'nCount_RNA')+ggtitle(i)+
          scale_y_log10()+scale_fill_manual(values = col_vector_2))
}
dev.off()

for(i in samples){
  Idents(seu[[i]]) <- 'seurat_clusters'
}

pdf(paste0(folder_results, "/3_filtered_UMAP_cluster.pdf"), height = 3.5, width = 4.5)
for(i in samples){
  print(UMAPPlot(object = seu[[i]], reduction = 'umap')+ggtitle(i))
}
dev.off()

pdf(paste0(folder_results, "/3_filtered_sctUMAP_cluster.pdf"), height = 3.5, width = 4.5)
for(i in samples){
  print(UMAPPlot(object = seu[[i]], reduction = 'SCTumap')+ggtitle(i))
}
dev.off()

##--------------------------------------------------------------------------------------------------##


##--------------------------------------------------------------------------------------------------##
saveRDS(seu, paste0(robject_folder, '0_filtered_normalised_individual_objects.RDS'))
##--------------------------------------------------------------------------------------------------##
