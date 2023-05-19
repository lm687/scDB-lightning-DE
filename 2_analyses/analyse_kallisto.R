##https://pachterlab.github.io/kallistobustools/tutorials/kb_getting_started/R/kb_intro_2_R/

##--------------------------------------------------------------------------------------------------##
rm(list = ls())
set.seed(1234)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

local <- F
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
folder_results <- "~/projects/lighting/3_results/" ## CCB cluster
robject_folder <- "~/projects/lighting/data/robjects/" ## CCB cluster
genome <- 'dmel649ChrimsonV2'
genome <- 'dmel649ChrimsonV2StringentFiltering' ## removing cells which have lower number of reads (but that still passed Charly's filter!)
input_objs <- paste0('KB-', genome)
# input_objs <- paste0('KB_onlyG1G2-', genome)

robject_folder <- paste0(robject_folder, input_objs, '/')
system(paste0("mkdir -p ", folder_results, input_objs))
system(paste0("mkdir -p ", robject_folder))

system(paste0("mkdir -p ", folder_results, input_objs, "/chrimsonpos/"))
system(paste0("mkdir -p ", folder_results, input_objs, "/chrimsonneg/"))
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
theme_set(theme_bw())
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## functions

source("helper_functions.R")

read_count_output <- function(dir, name) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- ".genes.txt"
  genes <- readLines(file(paste0(dir, "/", name, ge)))
  barcodes <- readLines(file(paste0(dir, "/", name, ".barcodes.txt")))
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

knee_plot <- function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  p <- ggplot(knee_plt, aes(total, rank)) +
    geom_line() +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2) +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    annotation_logticks() +
    labs(y = "Rank", x = "Total UMIs")
  return(p)
}

plot_gene_rank <- function(markers, n) {
  ## https://pachterlab.github.io/kallistobustools/tutorials/kb_building_atlas/R/kb_analysis_0_R/
  df_plot <- markers %>%
    group_by(cluster) %>%
    top_n(25, avg_log2FC) %>%
    mutate(rank = factor(row_number(desc(avg_log2FC))))
  ggplot(df_plot, aes(rank, avg_log2FC)) +
    geom_text(aes(label = gene), angle = -90, hjust = 1) +
    facet_wrap(~ cluster) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
}
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
samples <- apply(expand.grid(paste0('G', 1:4), paste0('_rep', 1:2)), 1, paste0, collapse='')
samples

if(grepl('onlyG1G2', input_objs)){
  samples <- samples[grepl('G1_', samples) | grepl('G2_', samples)]
}
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## read counts; check presence of chrimson
res_mat <- lapply(samples, function(i){
  read_count_output(paste0("~/projects/lighting/data/lightning/alignment/KallistoOUT/", gsub('StringentFiltering', '', genome), "/",
  i, "/counts_unfiltered"), name = "cells_x_genes")
})
names(res_mat) <- samples
chrimson_features_cells <- sapply(res_mat, function(i) table(i['FBto0000555',])) ## there is some chrimson i all samples!
sapply(chrimson_features_cells, function(i) 100*(1-(i[1]/sum(i)))) ## percentage of chrimson pos. it is very low...
## these results are the same as those with the dmel649Chrimson reference genome
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## filtering empty droplets

tot_counts <- lapply(res_mat, Matrix::colSums)
bc_rank <- lapply(res_mat, barcodeRanks, lower = 10)
names(tot_counts) <- samples
names(bc_rank) <- samples

## knee plots
options(repr.plot.width=9, repr.plot.height=6)
knee_plots <- lapply(samples, function(i) knee_plot(bc_rank[[i]])+ggtitle(i))
names(knee_plots) <- samples
knee_plots <- knee_plots[gtools::mixedorder(names(knee_plots))]
do.call('grid.arrange', knee_plots)
pdf(paste0(folder_results, input_objs, "/kneeplots.pdf"), height = 5.5, width = 3.5)
grid.arrange(grobs=knee_plots, nrow=4)
dev.off()
        
## filtering
res_mat <- lapply(samples, function(i) res_mat[[i]][, tot_counts[[i]] > metadata(bc_rank[[i]])$inflection])
res_mat <- lapply(res_mat, function(i) i[Matrix::rowSums(i) > 0,])
names(res_mat) <- samples
sapply(res_mat, dim)

##--------------------------------------------------------------------------------------------------##

## genes
name_conversion_file <- give_name_conversion_file()
name_conversion_file <- rbind(name_conversion_file, c('FBto0000555', 'Chrimson'))
names(name_conversion_file) <- c('gene', 'gene_symbol')
head(name_conversion_file)

for(i in samples){
  rownames(res_mat[[i]]) <- name_conversion_file$gene_symbol[match(rownames(res_mat[[i]]), name_conversion_file$gene)]
}
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## independent objects to be created after
seu <- lapply(res_mat, function(i) CreateSeuratObject(i, min.cells = 3, min.features = 200))

## remove cells if there is any further filtering that we want to do
if(genome == 'dmel649ChrimsonV2StringentFiltering'){
  warning('Filtering cells with further filters')
  Fullobj <- readRDS("~/projects/lighting/data/robjects/KB-dmel649ChrimsonV2/combined_dataset.RDS")
  Fullobj$origin <- gsub("^.*_","", colnames(Fullobj@assays$RNA@counts))
  StringentFilteringobj <- readRDS("~/projects/DrosophilaBrainAnnotation/out/robjects/annotation_Chrimson/KB-dmel649ChrimsonV2/chrimson_dataset.RDS")
  table(Fullobj$origin)
  table(StringentFilteringobj$origin)
  plot(as.vector(table(Fullobj$origin)), as.vector(table(StringentFilteringobj$origin)));abline(coef = c(0,1))

  dimnames(StringentFilteringobj)[[2]]
  StringentFilteringobj$origin <- samples[as.numeric(StringentFilteringobj$origin)]
  text(sapply(seu[samples], function(i) length(dimnames(i)[[2]])),
       as.vector(table(StringentFilteringobj$origin)[samples]), label=samples);abline(coef = c(0,1))
  ## remove the cells
  seu_nofilt <- seu
  for(i in samples){
    stopifnot(all(gsub("_.*", "", dimnames(StringentFilteringobj)[[2]][StringentFilteringobj$origin == i]) %in% dimnames(seu[[i]])[[2]]))
    seu[[i]] <- subset(seu[[i]], cells=gsub("_.*", "", dimnames(StringentFilteringobj)[[2]][StringentFilteringobj$origin == i]))
  }
  sapply(samples, function(i) c(nofilt=dim(seu_nofilt[[i]])[2], filt=dim(seu[[i]])[2]))
  # >   sapply(samples, function(i) c(nofilt=dim(seu_nofilt[[i]])[2], filt=dim(seu[[i]])[2]))
  # G1_rep1 G2_rep1 G3_rep1 G4_rep1 G1_rep2 G2_rep2 G3_rep2 G4_rep2
  # nofilt   13535   13649   11099    2205   44191   40246   16526    9843
  # filt      4323    3448    3943    1808    6942    6816    4993    3473
  rm(seu_nofilt)
}

for(i in samples){
  seu[[i]][["percent.mt"]] <- PercentageFeatureSet(seu[[i]], pattern = "^mt:")
}

ggplot(rownames_as_col(melt((sapply(samples, function(i) mean(seu[[i]][["percent.mt"]][,1], na.rm=T))))),
       aes(x=name, y=value, group=1))+geom_point()+geom_line()+
  theme(axis.text.x = element_text(angle = 30, hjust=1))+labs(x='Sample', y='Percentage of mt')
ggsave(paste0(folder_results, input_objs, "/percent_mt.pdf"), height = 3.5, width = 3.5)

options(repr.plot.width=12, repr.plot.height=6)

vlplots <- lapply(seu, function(i) VlnPlot(i, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1))
do.call('grid.arrange', vlplots, n)

ggplot(melt(lapply(seu, function(i) i[[c("nFeature_RNA", "nCount_RNA", "percent.mt")]])), aes(x=L1, y=value))+
  geom_violin()+facet_wrap(.~variable, scales='free', nrow=3)+scale_y_log10()+
  theme(axis.text.x = element_text(angle = 30, hjust=1))+labs(x='Sample')
ggsave(paste0(folder_results, input_objs, "/violinplots.pdf"), height = 5.5, width = 3.5)


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
samples <- gtools::mixedsort(samples)

try(dev.off())
pdf(paste0(folder_results, input_objs, "/UMI_ngenes.pdf"), height = 7.5, width = 7.5)
grid.arrange(grobs=plots_UMI_genesdetected, ncol=2)
dev.off()

plots_UMI_mt <- lapply(samples, function(i) ggplot(seu[[i]]@meta.data, aes(nCount_RNA, percent.mt)) +
  geom_pointdensity() +
  scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
  labs(x = "Total UMI counts", y = "Percentage mitochondrial")+ggtitle(i))
pdf(paste0(folder_results, input_objs, "/UMI_mt.pdf"), height = 7.5, width = 7.5)
grid.arrange(grobs=plots_UMI_mt, ncol=2)
dev.off()

plots_UMI_mt_v2 <- lapply(samples, function(i) ggplot(seu[[i]]@meta.data, aes(nCount_RNA, percent.mt)) +
                         geom_pointdensity() +
                         scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
                         #scale_x_log10()+scale_y_log10()+
                        lims(x=c(min(seu[[i]]@meta.data$nCount_RNA),500), y=c(0,10))+
                         labs(x = "Total UMI counts", y = "Percentage mitochondrial")+ggtitle(i))
pdf(paste0(folder_results, input_objs, "/UMI_mt_v2.pdf"), height = 7.5, width = 7.5)
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

seu <- lapply(seu, FindVariableFeatures, nfeatures = 3000)
top10 <-  lapply(seu, function(i) head(VariableFeatures(i), n=10))
plot1 <-  lapply(seu, VariableFeaturePlot, log = FALSE)
HVG_plot <- lapply(samples, function(i) LabelPoints(plot1[[i]], points = top10[[i]], repel = TRUE))

do.call('grid.arrange', HVG_plot)

ggplot(melt(!apply(firstcol_to_rownames(dcast(melt(top10), value~L1)), 2, is.na)) %>% dplyr::filter(value),
       aes(y=Var1, x=Var2, fill=value))+geom_tile()+labs(fill='Is Top10 HVG')+
  theme(axis.text.x = element_text(angle = 30, hjust=1))+
  scale_fill_manual(values = c("darksalmon"))+labs(x='Sample', y='HVG')+theme(legend.position = "bottom")
# scale_fill_jcolors(palette = 'pal5')
  # scale_fill_manual(values = c("darksalmon", "darkgoldenrod1"))
ggsave(paste0(folder_results, input_objs, "/HVG_top10_tile.pdf"), height = 6.5, width = 4.5)

##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## PCA

plot_elbow_vlines <- function(seurat_obj){
  pca_obj <- Reductions(seu$G1_rep1, 'pca')
  ggplot(cbind.data.frame(x=1:length(pca_obj@stdev), y=pca_obj@stdev))+
    geom_segment(aes(x=x,xend=x,y=0,yend=y))+geom_point(aes(x=x, y=y))+labs(x='PC', y='Standard deviations')
}

seu <- lapply(seu, function(i) RunPCA(i, verbose = FALSE, npcs = 20)) # uses HVG by default
# elbows <- lapply(samples, function(i) ElbowPlot(seu[[i]], type='h', ndims = 20)+ggtitle(i))
elbows <- lapply(samples, function(i) plot_elbow_vlines(seu[[i]])+ggtitle(i))
try(dev.off())
pdf(paste0(folder_results, input_objs, "/pca_elbow.pdf"), height = 4.5, width = 6.5)
grid.arrange(grobs=elbows, ncol=4)
dev.off()


pcas <- lapply(samples, function(i) PCAPlot(seu[[i]])+ggtitle(i))
pdf(paste0(folder_results, input_objs, "/pcas.pdf"), height = 7.5, width = 15)
grid.arrange(grobs=pcas, ncol=4)
dev.off()
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## find neighbours/clustering
seu <- lapply(seu, FindNeighbors, dims = 1:10)
seu <- lapply(seu, FindClusters)
pcas <- lapply(samples, function(i) PCAPlot(seu[[i]])+ggtitle(i))
pdf(paste0(folder_results, input_objs, "/pcas_clusters.pdf"), height = 7.5, width = 15)
grid.arrange(grobs=pcas, ncol=4)
dev.off()
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## T-sne
seu <- lapply(seu, RunTSNE, dims = 1:10)

tsnes <- lapply(samples, function(i) TSNEPlot(seu[[i]])+ggtitle(i))
pdf(paste0(folder_results, input_objs, "/tsne_clusters.pdf"), height = 7.5, width = 15)
grid.arrange(grobs=tsnes, ncol=4)
dev.off()


## umap
seu <- lapply(seu, RunUMAP, dims = 1:10)

umaps <- lapply(samples, function(i) UMAPPlot(seu[[i]])+ggtitle(i))
pdf(paste0(folder_results, input_objs, "/umap_clusters.pdf"), height = 7.5, width = 15)
grid.arrange(grobs=umaps, ncol=4)
dev.off()
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## biomarkers
biomarkers_list <- readRDS("/project/sims-lab/lmorrill/robjects/lightning/biomarkers_list.RDS")
biomarkers_list_vec <- c(biomarkers_list$`gene name`, 'Vmat', 'DAT', 'Frq1', 'otp', 'Lmx1a', 'Tdc2', 'Tbh', 'Tbh', 'kek1', 'hth', 'mirr', 'vvl', 'ey', 'Lim3')


## { nSyb, MtnA, Frq1, ey, Dop1R2, ct, elav, ..., several others } look very interesting to find the common clusters between samples
system(paste0("mkdir -p ", folder_results, input_objs, "/biomarkers/"))
for(bm in biomarkers_list_vec){
  cat(bm, which(biomarkers_list_vec == bm), "/", length(biomarkers_list_vec), '\n')
  try({
  umaps_bm <- lapply(samples, function(i) FeaturePlot(seu[[i]], features = bm)+ggtitle(paste0(bm, ' ', i)))
  pdf(paste0(folder_results, input_objs, "/biomarkers/umap_", bm, ".pdf"), height = 7.5, width = 15)
  grid.arrange(grobs=umaps_bm, ncol=4, main=bm)
  dev.off()
  try(rm(umaps_bm))
  })
}

pdf(paste0(folder_results, input_objs, "/biomarkers/umap_", 'Chrimson', ".pdf"), height = 7.5, width = 15)
grid.arrange(grobs=lapply(samples, function(i) FeaturePlot(seu[[i]], features = 'Chrimson')+ggtitle(paste0('Chrimson', ' ', i))),
             ncol=4, main=bm)
dev.off()

##--------------------------------------------------------------------------------------------------##


##--------------------------------------------------------------------------------------------------##
## integrate the datasets
## https://satijalab.org/seurat/articles/integration_introduction.html
#' These methods first identify cross-dataset pairs of cells that are in a matched biological state ('anchors'),
#' can be used both to correct for technical differences between datasets (i.e. batch effect correction), and to perform comparative scRNA-seq analysis of across experimental conditions.

features_for_integration <- SelectIntegrationFeatures(object.list = seu)
anchors_integration <- FindIntegrationAnchors(object.list = seu, anchor.features = features_for_integration)
combined_dataset <- IntegrateData(anchorset = anchors_integration)

save.image(paste0(robject_folder, 'individual_and_combined_datasets.RData'))
# load(paste0(robject_folder, 'individual_and_combined_datasets.RData'))
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Analyse the integrated datase

combined_dataset@assays$integrated@data == combined_dataset@assays$integrated@counts
# combined_dataset <- NormalizeData(combined_dataset, normalization.method = "LogNormalize", scale.factor = 10000) # https://github.com/satijalab/seurat/issues/4228
# FindVariableFeatures(combined_dataset, selection.method = "vst", nfeatures = 2000) ## we already have them - it's done on 2000 HVG
combined_dataset <- ScaleData(combined_dataset)

combined_dataset <- RunPCA(combined_dataset)
FeaturePlot(combined_dataset, reduction = "pca", feature = "Vmat")
options(repr.plot.width=9, repr.plot.height=6)
ElbowPlot(combined_dataset)

combined_dataset <- RunUMAP(combined_dataset, dims = 1:10, verbose = FALSE)
combined_dataset <- RunTSNE(combined_dataset, dims = 1:10, verbose = FALSE)
DimPlot(combined_dataset, reduction = "umap")

# pbmc <- FindNeighbors(pbmc, dims = 1:10, k.param = 20)
# pbmc <- FindClusters(pbmc, resolution = 0.6)

##--------------------------------------------------------------------------------------------------##

## compare with the cellranger dataset
## read in the cellranger datasets
FEB2019_noTE <- readRDS("~/projects/lighting/data/robjects/CharlyCellRanger/FEB2019_noTE.RDS")
## this object <FEB2019_noTE> had the objects integrated in this order: G1_rep1, G1_rep2, G2_rep1, etc.
name_sample_it <- 1

compare_alignments <- cbind.data.frame(cellranger=c('1_1',  '1_2',  '1_3',  '1_4',  '1_5',  '1_6',  '1_7',  '1_8'), kallisto=samples,
                                       stim=c('G1', 'G2', 'G3', 'G4', 'G1', 'G2', 'G3', 'G4'),
                                       condition=c('lights ON', 'lights OFF', 'lights ON', 'lights OFF', 'lights ON', 'lights OFF', 'lights ON', 'lights OFF'))

table(gsub(".*-", "", colnames(FEB2019_noTE@assays$RNA))) ## for each of the samples
# 1_1  1_2  1_3  1_4  1_5  1_6  1_7  1_8 
# 9268 9802 8814 9765 8931 9393 2138 8180 
  
give_cor_plots_CR_KB <- function(name_sample_it, plot_correlation=T, plot_chrimson=F, return_df=F){
  ## matching barcodes
  .mtch_CR_KB <- match(paste0(colnames(seu[[compare_alignments$kallisto[name_sample_it]]]@assays$RNA), '-', compare_alignments$cellranger[name_sample_it]),
                       colnames(FEB2019_noTE@assays$RNA))
  length(.mtch_CR_KB) == ncol(FEB2019_noTE@assays$RNA)
  stopifnot(length(.mtch_CR_KB) == ncol(seu[[compare_alignments$kallisto[name_sample_it]]]@assays$RNA))
  .subset_CR <- FEB2019_noTE@assays$RNA[,remove_na(.mtch_CR_KB)]
  .subset_KB <- seu[[compare_alignments$kallisto[name_sample_it]]]@assays$RNA[,!is.na(.mtch_CR_KB)]
  dim(.subset_CR)[2] == dim(.subset_KB)[2]
  
  ## matching genes
  .mtch_CR_KB2 <- match(rownames(.subset_CR),
                        rownames(.subset_KB))
  length(.mtch_CR_KB2) == nrow(.subset_CR)
  stopifnot(length(.mtch_CR_KB2) == nrow(.subset_CR))
  .subset_CR <- .subset_CR[!is.na(.mtch_CR_KB2),]
  .subset_KB <- .subset_KB[remove_na(.mtch_CR_KB2),]
  stopifnot(all(dim(.subset_CR) == dim(.subset_KB)))
  
  subset_vec_CR_KB <- sample(1:prod(dim(.subset_CR)), size = 2000, replace = F)
  
  if(plot_correlation){
    df <- data.frame(CR=as.vector(.subset_CR)[subset_vec_CR_KB], KB=as.vector(.subset_KB)[subset_vec_CR_KB])
    a <- ggplot(df,
                aes(CR, KB)) +
      geom_abline(slope = 1, intercept = 0, lty='dashed')+
      geom_pointdensity() +
      scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
      labs(x = "Counts Cell Ranger", y = "Counts Kallisto")+ggtitle(compare_alignments$kallisto[name_sample_it])+
      guides(col='none')
  }else if(plot_chrimson){
    df <- data.frame(CR=as.vector(.subset_CR['Chrimson',]), KB=as.vector(.subset_KB['Chrimson',]))
    a <- ggplot(df,
                aes(CR, KB)) +
      geom_abline(slope = 1, intercept = 0, lty='dashed')+
      geom_pointdensity() +
      scale_color_scico(palette = "devon", direction = -1, end = 0.9) +
      labs(x = "Chrimson counts Cell Ranger", y = "Chrimson counts Kallisto")+ggtitle(compare_alignments$kallisto[name_sample_it])+
      guides(col='none')
  }
  
  if(return_df){
    return(df)
  }else{
    return(a)
  }
}
plots_CR_KB_cor <- lapply(1:8, give_cor_plots_CR_KB)
plots_CR_KB_cor_chrimson <- lapply(1:8, give_cor_plots_CR_KB, plot_correlation = F, plot_chrimson = T)
plots_CR_KB_cor_chrimson_df <- lapply(1:8, give_cor_plots_CR_KB, plot_correlation = F, plot_chrimson = T, return_df=T)

## compare the raw counts
pdf(paste0(folder_results, input_objs, "/comparison_correlation_KB_CR", '', ".pdf"), height = 5.5, width = 7.5)
grid.arrange(grobs=plots_CR_KB_cor, nrow=2)
dev.off()

## compare which cells are chrimson+, by matching the barcode
pdf(paste0(folder_results, input_objs, "/comparison_correlation_Chrimson_KB_CR", '', ".pdf"), height = 5.5, width = 7.5)
grid.arrange(grobs=plots_CR_KB_cor_chrimson, nrow=2)
dev.off()

## venn diagrams
chrimsonpos_comparison <- sapply(1:8, function(i) table(apply(plots_CR_KB_cor_chrimson_df[[i]] > 0, 1, paste0, collapse='-')))
names(chrimsonpos_comparison) <- compare_alignments$kallisto
ggplot(melt(chrimsonpos_comparison), aes(x=L1, y=value, fill=Var1))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values = c('#c7fcf2', '#fc8897', 'red', 'green'))+
  labs(fill='Chrimson status (CR vs KB)', x='Samples', y='Number of cells')+
  theme(axis.text.x = element_text(angle = 30, hjust=1))+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave(paste0(folder_results, input_objs, "/comparison_barplot_Chrimson_KB_CR", '', ".pdf"), height = 4.5, width = 5.5)

##--------------------------------------------------------------------------------------------------##
## Analysis of the integrated dataset (cont)

DimPlot(combined_dataset, reduction = "umap")
for(bm in biomarkers_list_vec){
  cat(bm, which(biomarkers_list_vec == bm), "/", length(biomarkers_list_vec), '\n')
  try({
    FeaturePlot(combined_dataset, features = bm)+ggtitle(paste0(bm, ' ', "(integrated)"))
    ggsave(paste0(folder_results, input_objs, "/biomarkers/umap_integrated_", bm, ".pdf"),
           height = 5, width = 6)
    try(rm(umaps_bm))
  })
}

## heatmap of biomarkers
# DefaultAssay(combined_dataset) <- 'RNA'
# combined_dataset <- ScaleData(combined_dataset)
## even with few genes it takes too long - I have not run it
# DoHeatmap(combined_dataset, features = biomarkers_list_vec[!(biomarkers_list_vec %in% c('CG8177', 'CG4797'))])

DefaultAssay(combined_dataset) <- 'integrate'

VlnPlot(combined_dataset, features = biomarkers_list_vec[!(biomarkers_list_vec %in% c('CG8177', 'CG4797'))])

##--------------------------------------------------------------------------------------------------##
## re-calculate groups for integrate dataset -- I don't know why there are clusters already, possibly from the individual datasets?
DefaultAssay(combined_dataset) ## integrated
DefaultAssay(combined_dataset) <- 'integrated'
ElbowPlot(combined_dataset)
combined_dataset <- FindNeighbors(combined_dataset, dims = 1:10)
combined_dataset <- FindClusters(combined_dataset, dims = 1:10)
combined_dataset <- FindClusters(combined_dataset, dims = 1:10, resolution = 0.1)
combined_dataset <- FindClusters(combined_dataset, dims = 1:10, resolution = 0.01)

saveRDS(combined_dataset, paste0(robject_folder, 'combined_dataset.RDS'))
# combined_dataset <- readRDS(paste0(robject_folder, 'combined_dataset.RDS'))

##--------------------------------------------------------------------------------------------------##

unique(combined_dataset$integrated_snn_res.0.8)
combined_dataset$seurat_clusters <- combined_dataset$integrated_snn_res.0.8
Idents(combined_dataset) <- combined_dataset$seurat_clusters
DefaultAssay(combined_dataset) <- "RNA"
combined_dataset_markers <- FindAllMarkers(combined_dataset, test.use = "wilcox", only.pos = TRUE, 
                               min.pct = 0.25, logfc.threshold = 0.25)
saveRDS(combined_dataset_markers, paste0(robject_folder, 'combined_dataset_markers.RDS'))
# combined_dataset_markers <- readRDS(paste0(robject_folder, 'combined_dataset_markers.RDS'))
combined_dataset_markers

plot_gene_rank(combined_dataset_markers, 25)
ggsave(paste0(folder_results, input_objs, "/clusters_plot_gene_rank.pdf"), height = 4.5, width = 8.5)

col_vector <- readRDS("~/small_practical_robjects/col_vector3.RDS")
umap_single_facet_with_topmarker(seurat_obj = combined_dataset,
                                 markers = combined_dataset_markers,
                                 seurat_name_clusters = 'seurat_clusters')
ggsave(paste0(folder_results, input_objs, "/clusters_umap_single_facet_with_topmarker_all.pdf"), height = 5.5, width = 5.5)

UMAPPlot(combined_dataset)
PCAPlot(combined_dataset)
TSNEPlot(combined_dataset, label=T)
splitDimRedPlot(combined_dataset, group.by = 'seurat_clusters', plot_only_true=T)
# splitDimRedPlot(combined_dataset, group.by = 'seurat_clusters', plot_only_true=T, dim_red = 'pca')
# splitDimRedPlot(combined_dataset, group.by = 'seurat_clusters', plot_only_true=T, dim_red = 'tsne')
## why is this so incredibly split?

##--------------------------------------------------------------------------------------------------##

Chrimsonpos <- combined_dataset@assays$RNA@counts['Chrimson',] > 0
Chrimsonposcells <- subset(combined_dataset, cells = which(Chrimsonpos))
Chrimsonnegcells <- subset(combined_dataset, cells = which(!Chrimsonpos))
Chrimsonposcells ## 1506 samples
Chrimsonnegcells ## 81078 samples

#--------------------------------

Chrimsonposcells <- FindVariableFeatures(Chrimsonposcells)
Chrimsonposcells <- RunPCA(Chrimsonposcells, npcs = 60)
ElbowPlot(Chrimsonposcells)
Chrimsonposcells <- add_dimred(Chrimsonposcells, nPCs = 10, resolution_clusters = 1)
UMAPPlot(Chrimsonposcells)
ggsave(paste0(folder_results, input_objs, "/chrimsonpos/clusters_umap.pdf"), height = 4.5, width = 4.5)
FeaturePlot(Chrimsonposcells, features = 'DAT')
ggsave(paste0(folder_results, input_objs, "/chrimsonpos/umap_DAT.pdf"), height = 4.5, width = 4.5)
VlnPlot(Chrimsonposcells, features = 'DAT')
ggsave(paste0(folder_results, input_objs, "/chrimsonpos/VlnPlot_DAT.pdf"), height = 4.5, width = 4.5)

#--------------------------------

Chrimsonnegcells <- FindVariableFeatures(Chrimsonnegcells) ## we already have only 2000 HVGs in the integrated dataset
Chrimsonnegcells <- RunPCA(Chrimsonnegcells, npcs = 60)
ElbowPlot(Chrimsonnegcells)
Chrimsonnegcells <- add_dimred(Chrimsonnegcells, nPCs = 10, resolution_clusters = 1)
UMAPPlot(Chrimsonnegcells)
ggsave(paste0(folder_results, input_objs, "/chrimsonneg/clusters_umap.pdf"), height = 4.5, width = 4.5)
FeaturePlot(Chrimsonnegcells, features = 'DAT')
ggsave(paste0(folder_results, input_objs, "/chrimsonneg/umap_DAT.pdf"), height = 4.5, width = 4.5)
VlnPlot(Chrimsonnegcells, features = 'DAT')
ggsave(paste0(folder_results, input_objs, "/chrimsonneg/VlnPlot_DAT.pdf"), height = 4.5, width = 4.5)

##--------------------------------------------------------------------------------------------------##
## Plot number of cells of each condition in each cluster
Chrimsonposcells@meta.data$light <- (Chrimsonposcells@meta.data$stim %in% c('G1', 'G3'))
Chrimsonnegcells@meta.data$light <- (Chrimsonnegcells@meta.data$stim %in% c('G1', 'G3'))
combined_dataset@meta.data$light <- (combined_dataset@meta.data$stim %in% c('G1', 'G3'))

# add_number_cells_condition_per_cluster <- function(dataset, name_condition, cols_clusters=NULL){
#   stopifnot( length(unique(dataset@meta.data[,name_condition]) ) == 2) ## only two conditions
#   first_cond = unique(dataset@meta.data[,name_condition])[1]
#   if(is.null(cols_clusters)){
#     ## add all the cluster columns
#     cols_clusters <- python_like_select(colnames(dataset@meta.data), 'snn_res')
#   }
#   cat('Number of cluster columns:', length(cols_clusters), '\n')
#   for(i in cols_clusters){
#     cat('Cluster column:',i, '\n')
#     ## for each cluster, compute the fraction of cells in the first condition
#     clusters_list <- sort(unique(dataset@meta.data[,i]))
#     clusters_fracs <- sapply(clusters_list, function(cluster_idx){
#       .x <- subset(dataset, cells = which(dataset@meta.data[,i] == cluster_idx))
#       sum(.x@meta.data[,name_condition] == first_cond)/nrow(.x@meta.data)
#     })
#     dataset@meta.data[,paste0(i, 'FracCond1')] <- clusters_fracs[match(dataset@meta.data[,i], names(clusters_fracs))]
#   }
#   return(dataset)
# }
# add_number_cells_condition_per_cluster(Chrimsonnegcells, name_condition = 'light', cols_clusters='integrated_snn_res.0.01')

FeaturePlot(Chrimsonnegcells)

## Differential expression

combined_dataset$stim <- compare_alignments$stim[match(gsub(".*_", "1_", names(combined_dataset$orig.ident)), compare_alignments$cellranger)]
Chrimsonnegcells$stim <- compare_alignments$stim[match(gsub(".*_", "1_", names(Chrimsonnegcells$orig.ident)), compare_alignments$cellranger)]
Chrimsonposcells$stim <- compare_alignments$stim[match(gsub(".*_", "1_", names(Chrimsonposcells$orig.ident)), compare_alignments$cellranger)]

DefaultAssay(Chrimsonnegcells) <- 'RNA'
DefaultAssay(Chrimsonposcells) <- 'RNA'
DefaultAssay(combined_dataset) <- 'RNA'
give_DE_analysis_lightON_lightOFF(dataset = Chrimsonnegcells, dataset_name = 'chrimsonneg')
give_DE_analysis_lightON_lightOFF(dataset = Chrimsonposcells, dataset_name = 'chrimsonpos')
give_DE_analysis_lightON_lightOFF(dataset = combined_dataset, dataset_name = 'allcells', add_name='_res0p01', name_clusters='integrated_snn_res.0.01')
give_DE_analysis_lightON_lightOFF(dataset = combined_dataset, dataset_name = 'allcells', add_name='_res0p1', name_clusters='integrated_snn_res.0.1')
give_DE_analysis_lightON_lightOFF(dataset = Chrimsonnegcells, dataset_name = 'chrimsonneg', add_name='_res0p01', name_clusters='integrated_snn_res.0.01')
give_DE_analysis_lightON_lightOFF(dataset = Chrimsonposcells, dataset_name = 'chrimsonpos', add_name='_res0p01', name_clusters='integrated_snn_res.0.01')

names_datasets_DE <- c('chrimsonpos', 'chrimsonneg')
DE_lights_per_cluster <- lapply(names_datasets_DE, function(i){
  cat(i, '\n')
  if(local){
    x <- readRDS(paste0("/Users/lenamorrill/Documents/projects/lightning/github-repo-lightning-DE/3_results_local/objects/", input_objs, "/", i, "_DEgenes_per_cluster.RDS"))
  }else{
    fle <- paste0("~/projects/lighting/data/robjects/", input_objs, "/", i, "_DEgenes_per_cluster.RDS")
    cat(fle, '\n')
    x <- readRDS(fle)
  }
  add_genenames(x)
})
names(DE_lights_per_cluster) <- names_datasets_DE
DE_lights_per_cluster_molten <- (reshape2::melt(DE_lights_per_cluster, id.vars=c('gene', 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj')))
DE_lights_per_cluster_molten_top_common <- DE_lights_per_cluster_molten[DE_lights_per_cluster_molten$gene %in% names(tail(sort(table(DE_lights_per_cluster_molten$gene)), n=20)),]
DE_lights_per_cluster_molten_top_common$L1L2 <- apply(DE_lights_per_cluster_molten_top_common[,c('L1', 'L2')],1, paste0, collapse='-')
DE_lights_per_cluster_molten_top_common$signif <- DE_lights_per_cluster_molten_top_common$p_val_adj <= 0.05

ggplot(DE_lights_per_cluster_molten_top_common,
       aes(y=factor(gene, levels=DE_lights_per_cluster_molten_top_common %>% 
                      group_by(gene) %>% 
                      dplyr::summarise(median(avg_log2FC)) %>% 
                      plyr::arrange(`median(avg_log2FC)`) %>% dplyr::select(gene) %>% unlist),
           x=avg_log2FC, fill=L1, alpha=0.2), lty=3)+geom_density_ridges()+theme_bw()+geom_vline(xintercept = 0)+
  scale_fill_jcolors(palette = "pal6")+labs(x='Average log-fold change across cluster', y='Gene', fill='')+
  theme(legend.position = "bottom")+guides(alpha="none")+
  facet_wrap(.~signif, ncol=2)+ggtitle('LFC split by statistical differential abundance')
ggsave(paste0(folder_results, input_objs, "/lightONvsOFF_summarygenes_joy.pdf"), height = 5, width = 7)

DE_lights_per_cluster_molten_top_common_levels <- DE_lights_per_cluster_molten_top_common %>% group_by(gene) %>% summarise(m=mean(avg_log2FC)) %>% as.vector
DE_lights_per_cluster_molten_top_common$gene <- factor(DE_lights_per_cluster_molten_top_common$gene, 
                                                       levels=DE_lights_per_cluster_molten_top_common_levels$gene[order(DE_lights_per_cluster_molten_top_common_levels$m)])
DE_lights_per_cluster_molten_top_common$L1L2 <- factor(DE_lights_per_cluster_molten_top_common$L1L2, 
                                                       levels=rev(gtools::mixedsort(unique(DE_lights_per_cluster_molten_top_common$L1L2))))
ggplot(DE_lights_per_cluster_molten_top_common, aes(y=gene, x=L1L2, fill=avg_log2FC))+
  geom_tile()+
  geom_point(data = DE_lights_per_cluster_molten_top_common %>% filter(signif), aes(y=gene, x=L1L2), col='black')+
  scale_fill_viridis_b()+
  theme(axis.text.x = element_text(angle = 30, hjust=1))+
  labs(x='Cluster', y='Gene')+
  facet_grid (.~ L1, scales = "free_x", space = "free_x")
# facet_wrap(.~interaction(L1), scales='free_x')+
ggsave(paste0(folder_results, input_objs, "/lightONvsOFF_summarygenes_tile.pdf"), height = 5, width = 10)

##--------------------------------------------------------------------------------------------------##

saveRDS(Chrimsonposcells, paste0(robject_folder, 'Chrimsonposcells.RDS'))
# Chrimsonposcells <- readRDS(paste0(robject_folder, 'Chrimsonposcells.RDS'))
saveRDS(Chrimsonnegcells, paste0(robject_folder, 'Chrimsonnegcells.RDS'))
# Chrimsonnegcells <- readRDS(paste0(robject_folder, 'Chrimsonnegcells.RDS'))

##--------------------------------------------------------------------------------------------------##
## do we need integration via anchors for combining the samples? do simple logNorm counts suffice?
combined_dataset$stim <- gsub(".*_", "G", names(combined_dataset$orig.ident))

infiles10X_merged <- merge(infiles10X[[1]],  y = infiles10X[-1], add.cell.ids = names_samples, project = 'GSE184507')

combined_dataset_merged <- merge(seu[[1]],  y = seu[-1], add.cell.ids = samples, project = 'KallistoChrimson')
combined_dataset_merged$stim <- sapply(names(combined_dataset_merged$orig.ident), function(i) paste0(strsplit(i, '_')[[1]][1:2], collapse = '_'))

DefaultAssay(combined_dataset_merged) <- 'RNA'
combined_dataset_merged <- FindVariableFeatures(combined_dataset_merged)
combined_dataset_merged <- ScaleData(combined_dataset_merged)
combined_dataset_merged <- RunPCA(combined_dataset_merged, npcs=20)
plot_elbow_vlines(combined_dataset_merged)
combined_dataset_merged <- RunUMAP(combined_dataset_merged, reduction = "pca", dims = 1:10, umap.method='uwot')

DefaultAssay(combined_dataset) <- 'integrated'
umap_int <- UMAPPlot(combined_dataset, group.by='stim') ## this is the dataset with integration
DefaultAssay(combined_dataset) <- 'RNA'

umap_not_int <- UMAPPlot(combined_dataset_merged, group.by='stim') ## this is the dataset without integration

try(dev.off())
pdf(paste0(folder_results, input_objs, "/comparison_integrated_notintegrated.pdf"), height = 5, width = 9)
cowplot::plot_grid(umap_int+ggtitle('Integrated'), umap_not_int+ggtitle('Not integrated'))
dev.off()
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Hr38: values per group
DEGs_df <- cbind.data.frame(seurat_clusters=combined_dataset[[]]$seurat_clusters,
                            integrated_snn_res.0.8=combined_dataset[[]]$integrated_snn_res.0.8,
                 stim=combined_dataset[[]]$stim,
                 Hr38=as.vector(combined_dataset@assays$RNA['Hr38',]),
                 sr=as.vector(combined_dataset@assays$RNA['sr',]),
                 CG14186=as.vector(combined_dataset@assays$RNA['CG14186',]))
sort_first_by_second <- function(i){
  as.vector(i[order(i[,2]),1])
}


ggplot((DEGs_df), aes(x=factor(seurat_clusters, levels=sort_first_by_second(DEGs_df %>% 
                                                             group_by(seurat_clusters) %>% summarize(mean(Hr38)) %>% 
                                                             as.data.frame)),
                                 y=Hr38+1, col=factor(stim, levels=paste0('G', c(1,3,2,4))), shape=(Hr38 == 0)))+
  geom_violin(position = position_dodge(width = 0.3))+
  scale_y_continuous(trans = "log2")+
  scale_color_manual(values=c('#f8938d', '#f49f57', '#328b5b', '#5c836d'))

basic_plot_DEGs_exprs <- ggplot(melt(DEGs_df[,!(colnames(DEGs_df) == "integrated_snn_res.0.8")],
                                     id.vars=c('seurat_clusters', 'stim')),
       aes(x=seurat_clusters,
                      y=value+1, fill=factor(stim, levels=paste0('G', c(1,3,2,4))), shape=(value == 0)))+
  scale_y_continuous(trans = "log2")+
  scale_fill_manual(values=c('#f8938d', '#f49f57', '#328b5b', '#34c7b6'))+
  facet_wrap(.~variable, ncol=1)+
  theme(legend.position = "bottom")+
  labs(fill='Group', x='Seurat cluster', y='Expression (counts)', shape='Zero counts')

basic_plot_DEGs_exprs+geom_violin(position = position_dodge(width = 0.3))
basic_plot_DEGs_exprs+geom_boxplot(position = position_dodge(width = 0.3))
ggsave(paste0(folder_results, input_objs, "/expression_ARGs_boxplot.pdf"), height = 6, width = 7.5)

basic_plot_DEGs_exprs_2 <- ggplot(melt(DEGs_df[,!(colnames(DEGs_df) == "seurat_clusters")],
                                     id.vars=c('integrated_snn_res.0.8', 'stim')),
                                aes(x=`integrated_snn_res.0.8`,
                                    y=value+1, fill=factor(stim, levels=paste0('G', c(1,3,2,4))), shape=(value == 0)))+
  scale_y_continuous(trans = "log2")+
  scale_fill_manual(values=c('#f8938d', '#f49f57', '#328b5b', '#34c7b6'))+
  facet_wrap(.~variable, ncol=1)+
  theme(legend.position = "bottom")+
  labs(fill='Group', x='Seurat cluster', y='Expression (counts)', shape='Zero counts')
basic_plot_DEGs_exprs_2+geom_violin(position = position_dodge(width = 0.8))


##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## annotate largest clusters
UMAPPlot(combined_dataset)
all(combined_dataset[[]]$seurat_clusters == combined_dataset[[]]$integrated_snn_res.0.01)
DefaultAssay(combined_dataset) <- 'RNA'
Idents(combined_dataset) <- (combined_dataset$integrated_snn_res.0.01)
markers0p01 <- Seurat::FindAllMarkers(combined_dataset, group.by='integrated_snn_res.0.01')
Idents(combined_dataset) <- (combined_dataset$integrated_snn_res.0.1)
markers0p1 <- Seurat::FindAllMarkers(combined_dataset, group.by='integrated_snn_res.0.1')

saveRDS(markers0p01, paste0(robject_folder, 'markers0p01.RDS'))
saveRDS(markers0p1, paste0(robject_folder, 'markers0p1.RDS'))
# markers0p01 <- readRDS(paste0(robject_folder, 'markers0p01.RDS'))
# markers0p1 <- readRDS(paste0(robject_folder, 'markers0p1.RDS'))


gen_marker_table <- function(x, markers0p01, n_arg=4){
  markers0p01[markers0p01$cluster == x, ] %>%
    head(n=n_arg)
}

TSNEPlot(combined_dataset)
top_markers_markers0p01 <- map_dfr(levels(markers0p01$cluster), gen_marker_table, markers0p01)
top_markers_markers0p01_top1 <- map_dfr(levels(markers0p01$cluster), gen_marker_table, markers0p01, n=1)
top_markers_markers0p1 <- map_dfr(levels(markers0p1$cluster), gen_marker_table, markers0p1)
top_markers_markers0p1_top1 <- map_dfr(levels(markers0p1$cluster), gen_marker_table, markers0p1, n=1)
top_markers_markers0p1_top20 <- map_dfr(levels(markers0p1$cluster), gen_marker_table, markers0p1, n=20)
# featureplots <- lapply(levels(top_markers_markers0p01$cluster), function(cl){ FeaturePlot(object = combined_dataset, 
#             features = c(top_markers_markers0p01[top_markers_markers0p01$cluster == cl, "gene"]), 
#             cols = c("grey", "#75c08c"), 
#             reduction = "umap", ncol=4) })
# grid.arrange(grobs=featureplots, nrow=4)
# do.call('grid.arrange', featureplots)

FeaturePlot(object = combined_dataset, 
           features = as.vector(sapply(levels(top_markers_markers0p01$cluster), function(cl){
             c(top_markers_markers0p01[top_markers_markers0p01$cluster == cl, "gene"])
             })), 
           cols = c("#eaa661", "#75c08c"), 
           reduction = "umap", ncol=4)

cowplot::plot_grid(UMAPPlot(combined_dataset),
                   FeaturePlot(object = combined_dataset, 
            features = as.vector(sapply(levels(top_markers_markers0p01$cluster), function(cl){
              c(top_markers_markers0p01[top_markers_markers0p01$cluster == cl, "gene"])
            }))[1], 
            cols = c("#eaa661", "#75c08c"), 
            reduction = "umap"))

top_markers_markers0p01_top1
FeaturePlot(object = combined_dataset, 
            features = 'rut',
            cols = c("#eaa661", "#75c08c"), 
            reduction = "umap")

rndm_subset_combined_dataset <- subset(combined_dataset, cells=sample(dimnames(combined_dataset)[[2]], size = 1000))
DoHeatmap(rndm_subset_combined_dataset, features = top_markers_markers0p01$gene)
DoHeatmap(rndm_subset_combined_dataset, features = top_markers_markers0p1_top20$gene)

UMAPPlot(combined_dataset, group.by='integrated_snn_res.0.1')
FeaturePlot(object = combined_dataset, 
            features = as.vector(sapply(levels(top_markers_markers0p1_top1$cluster), function(cl){
              c(top_markers_markers0p1_top1[top_markers_markers0p1_top1$cluster == cl, "gene"])
            })), 
            cols = c("#eaa661", "#75c08c"), 
            reduction = "umap", ncol=4)

cowplot::plot_grid(UMAPPlot(combined_dataset, group.by='integrated_snn_res.0.1'),
                   FeaturePlot(object = combined_dataset, 
                               features = as.vector(sapply(levels(top_markers_markers0p1_top1$cluster), function(cl){
                                 c(top_markers_markers0p1_top1[top_markers_markers0p1_top1$cluster == cl, "gene"])
                               }))[1], 
                               cols = c("#eaa661", "#75c08c"), 
                               reduction = "umap", ncol=4)
)

runif(10)

## annotation seems super difficult in the kallisto dataset - could it be a problem of what data is used to create
## the groups, or the data that is plotted?
FEB2019_noTE$seurat_clusters

##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
library(scran)
##' I think I should be using log, but not lognorm, data, for this :
##' For getDenoisedPCs, a numeric matrix of log-expression values, where rows
##' are genes and columns are cells. Alternatively, a SummarizedExperiment
##' object containing such a matrix. 
##' Moreover, are these data scaled?


combined_dataset_sce <- SingleCellExperiment(assays=list(counts=combined_dataset@assays$RNA@counts))
combined_dataset_sce <- scuttle::logNormCounts(combined_dataset_sce)
modelGeneVarres <- scran::modelGeneVar(x = combined_dataset_sce@assays@data$logcounts, block=combined_dataset$stim)
getTopHVGsres <- scran::getTopHVGs(combined_dataset_sce)
denoisedPCA <- scran::denoisePCA(x = combined_dataset_sce, technical=modelGeneVarres$tech, subset.row=getTopHVGsres)

colData(denoisedPCA)[,'Group'] <- gsub(".*_", "", rownames(colData(denoisedPCA)))

Reductions(denoisedPCA)
plotReducedDim(denoisedPCA)
plotPCA(denoisedPCA)
plotPCA(denoisedPCA, colour_by = 'Group')

## with decontex, find if any cells are contamination
library(celda)
decontx_res <- celda::decontX(denoisedPCA)
decontx_res
umap_decontX <- reducedDim(decontx_res, "decontX_UMAP")
plotDimReduceCluster(x = decontx_res$decontX_clusters, dim1 = umap_decontX[, 1], dim2 = umap_decontX[, 2])
plotDecontXContamination(decontx_res) ## really high contamination!
ggsave(paste0(folder_results, input_objs, "/decontX_contamination_UMAP.pdf"), height = 5, width = 6)

FeaturePlot(decontx_res, 'decontX_UMAP')

## detect doublets
library(scDblFinder)
#' Warning message:
#'   In scDblFinder(denoisedPCA) :
#'   You are trying to run scDblFinder on a very large number of cells. 
#'   If these are from different captures, please specify this using the `samples` argument.TRUE
sce_scdblfinder <- scDblFinder(denoisedPCA, samples = colData(denoisedPCA)$Group)
sce_scdblfinder

denoisedPCA <- scater::runUMAP(denoisedPCA)

ggplot(data.frame(reducedDim(sce_scdblfinder, "PCA")[,1:2], col=colData(sce_scdblfinder)$scDblFinder.class),
       aes(x=PC1, y=PC2, col=col))+
  geom_point()

table(colData(sce_scdblfinder)$scDblFinder.class)

ggplot(data.frame(reducedDim(denoisedPCA , "UMAP")[,1:2], col=colData(sce_scdblfinder)$scDblFinder.class),
       aes(x=X1, y=X2, col=col))+
  geom_point(size=0.2)+
  scale_color_manual(values = c("#fdb08f", "#8ee2c4"))+
  labs(x='UMAP 1', y='UMAP 2')
ggsave(paste0(folder_results, input_objs, "/doublets_UMAP.pdf"), height = 5, width = 6)

denoisedPCA
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## Cluster-unaware DE with Milo
library(miloDE)

combined_dataset@meta.data[,'stim_light'] <- combined_dataset[[]]$stim %in% c('G1', 'G3')
table(combined_dataset@meta.data[,'stim_light'] )

combined_dataset_sce <- scater::runPCA(combined_dataset_sce) ## PCA needs to be run for Milo
combined_dataset_sce <- scater::runUMAP(combined_dataset_sce) ## there is one cluster of only non-light cells
colData(combined_dataset_sce)[,'stim_light'] <- combined_dataset[[]][,'stim_light']
colData(combined_dataset_sce)[,'sample_id'] <- gsub(".*_", "", colnames(combined_dataset_sce))
head(colData(combined_dataset_sce))

ggplot(data.frame(reducedDim(combined_dataset_sce, "PCA")[,1:2], col=colData(combined_dataset_sce)$stim_light),
       aes(x=PC1, y=PC2, col=col))+
  geom_point()
ggplot(data.frame(reducedDim(combined_dataset_sce , "UMAP")[,1:2], col=colData(combined_dataset_sce)$stim_light),
       aes(x=X1, y=X2, col=col))+
  geom_point(size=0.2)+
  scale_color_manual(values = c("#fdb08f", "#8ee2c4"))+
  labs(x='UMAP 1', y='UMAP 2')


combined_dataset_sce = assign_neighbourhoods(combined_dataset_sce, k = 20, order = 2, 
                                        filtering = TRUE, reducedDim_name = "PCA")
de_stat = de_test_neighbourhoods(combined_dataset_sce , sample_id = "sample_id", 
                                 design = ~stim_light, covariates='stim_light') ## there is a covariates argument
de_stat_Hr38 <- de_stat %>% filter(gene=='Hr38')
apply(combined_dataset_sce@nhoods, 1, function(i) which(i == 1)) ## one cell can belong to several neighbourhoods
# colData(combined_dataset_sce)$nhood_id <- 
# nhoods_Hr38_DE <- combined_dataset_sce@nhoods[,remove_na(de_stat_Hr38$Nhood[de_stat_Hr38$pval_corrected_across_genes < 0.05])]
## select neighbourhoods with DE
# nhoods_Hr38_DE <- combined_dataset_sce@nhoods[,remove_na(de_stat_Hr38$Nhood[de_stat_Hr38$pval_corrected_across_nhoods < 0.05])]
# nhoods_Hr38_DE <- nhoods_Hr38_DE[rowSums(nhoods_Hr38_DE)>0,]
# dim(nhoods_Hr38_DE)
# de_stat_Hr38$Nhood[de_stat_Hr38$pval_corrected_across_genes < 0.05]
# colData(combined_dataset_sce)[,'Hr38_DE'] <- 

saveRDS(combined_dataset_sce, paste0(robject_folder, 'combined_dataset_sce.RDS'))
saveRDS(de_stat, paste0(robject_folder, 'de_stat.RDS'))


p1 = plot_milo_by_single_metric(sce_milo, stat_de_magnitude, colour_by = "n_DE_genes" , 
                                layout = "UMAP" , size_range = c(1.5,3) , edge_width = c(0.2,0.5)) + 
  scale_fill_viridis(name = "# DE genes")

de_stat

##--------------------------------------------------------------------------------------------------##
