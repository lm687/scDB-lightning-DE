##--------------------------------------------------------------------------------------------------##
rm(list = ls())
source("projects/lighting/2_analyses/helper_functions.R")
theme_set(theme_bw())

##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
library(SingleR)
library(clustree)
library(plyr)
library(ggalluvial)
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
robject_folder <- "~/projects/lighting/data/robjects/" ## CCB cluster
folder_results <- "~/projects/lighting/3_results/" ## CCB cluster
genome <- 'dmel649ChrimsonV2'
input_objs <- paste0('KB-', genome)
robject_folder <- paste0(robject_folder, input_objs, '/')

combined_dataset_light <- readRDS(paste0(robject_folder, 'combined_dataset.RDS'))
combined_dataset_lightCharly <- readRDS("~/projects/lighting/data/robjects/CharlyCellRanger/FEB2019_noTE.RDS")
aerstsseurat <- readRDS("~/projects/general/fly_midbrain_atlases/Atlas/AerstsSeurat.RDS")
aerstsseurat_consensus_single <- readRDS("~/projects/general/fly_midbrain_atlases/Atlas/AerstsSeuratConsensusCells.RDS")

##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## SingleR for cluster annotation
##' de.method="wilcox"
##'  we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise
##'  comparison between labels. This is slower but more appropriate for single-cell data
##'  compared to the default marker detection algorithm (which may fail for low-coverage
##'  data where the median is frequently zero).

## I haven't been able to finish running this without subsetting
# prediction_light_dataset_Aerts <- SingleR(test = GetAssayData(combined_dataset_light)[1:2000,1:2000],
#                                           ref = GetAssayData(aerstsseurat)[1:2000,1:2000],
#                                           labels = aerstsseurat@meta.data$annotation[1:2000])

## using a single cell from the Aerts clusters (expression computed using AverageExpression)
predict_with_Aerts_consensusCell <- function(test_assay){
  SingleR(test = test_assay,
          ref = aerstsseurat_consensus_single,
          labels = colnames(aerstsseurat_consensus_single))
  }
##--------------------------------------------------------------------------------------------------##

##--------------------------------------------------------------------------------------------------##
## predict

DefaultAssay(combined_dataset_light) <- 'RNA'
aerstsseurat ## 17473 features across 56902 samples within 1 assay 
combined_dataset_light # 14770 features across 82584 samples within 2 assays 
stopifnot(DefaultAssay(combined_dataset_light) == DefaultAssay(aerstsseurat))
## make sure that both are lognorm
prediction_light_dataset_Aerts_consensusCell <- predict_with_Aerts_consensusCell(GetAssayData(combined_dataset_light))
saveRDS(prediction_light_dataset_Aerts_consensusCell, paste0(robject_folder, '/annotation_cells/combined_dataset.RDS'))
# prediction_light_dataset_Aerts_consensusCell <- readRDS(paste0(robject_folder, '/annotation_cells/combined_dataset.RDS'))

# plotScoreHeatmap(prediction_light_dataset_Aerts)


## Charlie's dataset
DefaultAssay(combined_dataset_lightCharly) <- 'RNA'
combined_dataset_lightCharly@assays$RNA@data ## logNorm
stopifnot(DefaultAssay(combined_dataset_lightCharly) == DefaultAssay(aerstsseurat))

combined_dataset_lightCharly
## takes approx 5 min
prediction_lightCharly_dataset_Aerts_consensusCell <- predict_with_Aerts_consensusCell(GetAssayData(combined_dataset_lightCharly))
saveRDS(prediction_lightCharly_dataset_Aerts_consensusCell, paste0('~/projects/lighting/data/robjects/CharlyCellRanger/', '/annotation_cells/combined_dataset.RDS'))

##--------------------------------------------------------------------------------------------------##

combined_dataset_light$annotationsingleR <- prediction_light_dataset_Aerts$pruned.labels[match(dimnames(combined_dataset_light)[[2]], rownames(prediction_light_dataset_Aerts))]
combined_dataset_light$annotationsingleRConsensusCell <- prediction_light_dataset_Aerts_consensusCell$pruned.labels[match(dimnames(combined_dataset_light)[[2]], rownames(prediction_light_dataset_Aerts_consensusCell))]
table(is.na(combined_dataset_light$annotationsingleR))
UMAPPlot(combined_dataset_light, group.by='annotationsingleR')


sort(table(combined_dataset_light$annotationsingleR))
umap_plot_single_group(combined_dataset_light, 'annotationsingleR', 'L3')
umap_plot_single_group(combined_dataset_light, 'annotationsingleR', 'L3', split_umap=T, na_as_other = T)
umap_plot_single_group(combined_dataset_light, 'annotationsingleR', 'Tyraminergic', split_umap=T, na_as_other = T)
umap_plot_single_group(combined_dataset_light, 'annotationsingleR', 'Astrocyte-like', split_umap=T, na_as_other = T)
umap_plot_single_group(combined_dataset_light, 'annotationsingleR', 'Ensheathing_glia', split_umap=T, na_as_other = T)

pdf(paste0(folder_results, input_objs, "/annotation_cells/UMAP_annotation_annotationsingleRConsensusCell.pdf"), height = 3, width = 3)
for(class_it in sort(unique(combined_dataset_light$annotationsingleRConsensusCell))){
  cat(class_it, '...\n')
  print(umap_plot_single_group(combined_dataset_light, 'annotationsingleRConsensusCell', class_it,
                               split_umap=T, na_as_other = T, show_other=F, plot_legend=F)+labs(x='UMAP 1', y='UMAP 2'))
}
dev.off()

##--------------------------------------------------------------------------------------------------##
combined_dataset_lightCharly$annotationsingleRConsensusCell <- prediction_lightCharly_dataset_Aerts_consensusCell$pruned.labels[match(dimnames(combined_dataset_lightCharly)[[2]],
                                                                                                                                      rownames(prediction_lightCharly_dataset_Aerts_consensusCell))]

umap_plot_single_group(combined_dataset_lightCharly, 'annotationsingleRConsensusCell', 'Ensheathing_glia', split_umap=T, na_as_other = T)

pdf(paste0(folder_results, 'CharlyCellRanger/', "/annotation_cells/UMAP_annotation_annotationsingleRConsensusCell.pdf"), height = 3, width = 3)
for(class_it in gtools::mixsedsort(unique(combined_dataset_lightCharly$annotationsingleRConsensusCell))){
  cat(class_it, '...\n')
  print(umap_plot_single_group(combined_dataset_lightCharly, 'annotationsingleRConsensusCell', class_it,
                               split_umap=T, na_as_other = T, show_other=F, plot_legend=F)+labs(x='UMAP 1', y='UMAP 2'))
}
dev.off()

##--------------------------------------------------------------------------------------------------##
# 
# prediction_lightCharly_dataset_Aerts_consensusCell <- readRDS(paste0('~/projects/lighting/data/robjects/CharlyCellRanger/', '/annotation_cells/combined_dataset.RDS'))
# prediction_light_dataset_Aerts_consensusCell <- readRDS(paste0(robject_folder, '/annotation_cells/combined_dataset.RDS'))
# combined_dataset_light$annotationsingleRConsensusCell <- prediction_light_dataset_Aerts_consensusCell$pruned.labels[match(dimnames(combined_dataset_light)[[2]], rownames(prediction_light_dataset_Aerts_consensusCell))]
# combined_dataset_lightCharly$annotationsingleRConsensusCell <- prediction_lightCharly_dataset_Aerts_consensusCell$pruned.labels[match(dimnames(combined_dataset_lightCharly)[[2]],
#                                                                                                                                       rownames(prediction_lightCharly_dataset_Aerts_consensusCell))]

## compare cell annotations
# head(names(combined_dataset_lightCharly$orig.ident))
# unique(gsub(".*-", "", names(combined_dataset_lightCharly$orig.ident)))
# unique(gsub(".*_", "", names(combined_dataset_light$orig.ident)))

dimnames_combined_dataset_lightCharly <- gsub('-1', '', dimnames(combined_dataset_lightCharly)[[2]])
mcth_light_datasets <- match(dimnames_combined_dataset_lightCharly, dimnames(combined_dataset_light)[[2]])
length(mcth_light_datasets)
length(dimnames_combined_dataset_lightCharly)
length(dimnames(combined_dataset_light)[[2]])

combined_dataset_light_df <- cbind.data.frame(charly=combined_dataset_lightCharly@meta.data[,c('seurat_clusters', 'annotationsingleRConsensusCell')],
                 kallisto=combined_dataset_light@meta.data[mcth_light_datasets,c('seurat_clusters', 'annotationsingleRConsensusCell')])

combined_dataset_light_confounding <- table(combined_dataset_light_df$charly.annotationsingleRConsensusCell,
      combined_dataset_light_df$kallisto.annotationsingleRConsensusCell)

pheatmap::pheatmap(combined_dataset_light_confounding)

combined_dataset_light_confounding_melt <- (cbind(melt(combined_dataset_light_confounding), frac=melt(normalise_rw(combined_dataset_light_confounding))[,3]))

combined_dataset_light_confounding_melt$Var1 <- factor(as.character(combined_dataset_light_confounding_melt$Var1), levels = gtools::mixedsort(unique(c(as.character(combined_dataset_light_confounding_melt$Var1),
                                                                                                                                                       as.character(combined_dataset_light_confounding_melt$Var2)))))
combined_dataset_light_confounding_melt$Var2 <- factor(as.character(combined_dataset_light_confounding_melt$Var2), levels = gtools::mixedsort(unique(c(as.character(combined_dataset_light_confounding_melt$Var1),
                                                                                                                                                       as.character(combined_dataset_light_confounding_melt$Var2)))))

ggplot(combined_dataset_light_confounding_melt %>% dplyr::filter(frac>0), aes(x=Var1, y=Var2,size=frac, col=value, shape=(as.character(Var1)==as.character(Var2))))+
  geom_point()+scale_colour_viridis_c()+theme_bw()+
  theme(axis.text.x = element_text(angle = 30, hjust=1))+
  labs(shape='Match', col='Num observations')
ggsave(paste0(folder_results, input_objs, "/annotation_cells/comparison_annotation_annotationsingleRConsensusCell_CharlyCellRanger.pdf"),
       height = 13, width = 14.5)

##--------------------------------------------------------------------------------------------------##

color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
color200 <- sample(color, size = 200, replace = F)
## do our clusters correspond with a single cell type?

ggplot(combined_dataset_light@meta.data[,c('seurat_clusters', 'annotationsingleRConsensusCell')],
       aes(x=seurat_clusters, fill=annotationsingleRConsensusCell))+
  geom_bar()


ggplot(combined_dataset_light@meta.data[,c('integrated_snn_res.0.8', 'annotationsingleRConsensusCell')],
       aes(x=integrated_snn_res.0.8, fill=annotationsingleRConsensusCell))+
  geom_bar()+facet_wrap(.~integrated_snn_res.0.8, scale='free', nrow=1)+
  theme(legend.position = "bottom")+guides(fill=guide_legend(ncol=10))+
  scale_fill_manual(values = color200)
ggsave(paste0(folder_results, input_objs, "/annotation_cells/comparison_clusters_annotationsingleRConsensusCell.pdf"),
       height = 5, width = 14.5)

ggplot(combined_dataset_lightCharly@meta.data[,c('integrated_snn_res.4', 'annotationsingleRConsensusCell')],
       aes(x=integrated_snn_res.4, fill=annotationsingleRConsensusCell))+
  geom_bar()+facet_wrap(.~integrated_snn_res.4, scale='free', nrow=7)+
  theme(legend.position = "bottom")+guides(fill=guide_legend(ncol=10))+
  scale_fill_manual(values = color200)
ggsave(paste0(folder_results, 'CharlyCellRanger', "/annotation_cells/comparison_clusters_annotationsingleRConsensusCell.pdf"),
       height = 15, width = 14.5)

## give fraction of cells in my clusters annotated to the most common category
factors_first_col_by_second <- function(i){
  i[,1] <- factor(i[,1], levels=as.vector(i[order(i[,2], decreasing = T),1]))
  i
}

## fraction of cells that belong to the most common annotated group, sorted
give_frac_cells_cluster_in_largest_group <- function(DTASET){
  ggplot(factors_first_col_by_second(data.frame(DTASET@meta.data %>% dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(max_frac=max(table(annotationsingleRConsensusCell)/sum(table(annotationsingleRConsensusCell)))))),
    aes(x=seurat_clusters, y=max_frac))+
    geom_point()+geom_line(aes(group=1))+
    theme_bw()
}

give_frac_cells_cluster_in_largest_group(combined_dataset_light)
give_frac_cells_cluster_in_largest_group(combined_dataset_lightCharly)

## scatterplot of the fraction of cells in each cluster that belong to the most common annotated group, and the size of the cluster
give_frac_cells_cluster_in_largest_group_size_scatter <- function(DTASET){
  ggplot(data.frame(DTASET@meta.data %>% dplyr::group_by(seurat_clusters) %>%
             dplyr::summarise(max_frac=max(table(annotationsingleRConsensusCell)/sum(table(annotationsingleRConsensusCell))),
                              size_cluster=n())), aes(x=size_cluster, y=max_frac))+
  geom_point()
}

give_frac_cells_cluster_in_largest_group_size_scatter(combined_dataset_lightCharly) ## negative correlation: larger clusters are the dirtiest
give_frac_cells_cluster_in_largest_group_size_scatter(combined_dataset_light)

give_table_matches <- function(combined_dataset_light_df, col1, col2, prune=F, prune_val=0.6){
  .x <- table(combined_dataset_light_df[,col1],
                                              combined_dataset_light_df[,col2])
  rownames(.x) <- paste0('Cluster', rownames(.x))
  
  .x_melt <- (cbind(melt(.x), frac=melt(normalise_rw(.x))[,3]))
  .x_melt$Var1 <- as.character(.x_melt$Var1)
  .x_melt$Var2 <- as.character(.x_melt$Var2)
  
  if(prune){
    .x_melt <- .x_melt %>% dplyr::filter(frac>prune_val)
  }else{
    .x_melt <- .x_melt %>% dplyr::filter(frac>0)
  }
  
  
  .x_melt[.x_melt$Var1 == 'Cluster3',]
  sort(normalise_rw(.x['Cluster3',]), decreasing = T)[1:5]
  sort(.x['Cluster3',], decreasing = T)[1:5]
  
  
  # .x_melt$Var1 <- factor(as.character(.x_melt$Var1), levels=gtools::mixedsort(unique(as.character(.x_melt$Var1))))
  
  .x_melt$Var1 <- factor(.x_melt$Var1, levels=gtools::mixedsort(unique(.x_melt$Var1)))
  .x_melt$Var2 <- factor(.x_melt$Var2, levels=unique(c(sapply(as.character(levels(.x_melt$Var1)), function(i) names(which.max(.x[i,]))),
                                   colnames(.x))))
  
  ggplot(.x_melt, aes(x=Var1, y=Var2,size=frac, col=value, shape=(as.character(Var1)==as.character(Var2))))+
    geom_point()+scale_colour_viridis_c()+theme_bw()+
    theme(axis.text.x = element_text(angle = 30, hjust=1))+
    labs(shape='Match', col='Num observations')
   
}

give_table_matches(combined_dataset_light_df = combined_dataset_light[[]], col1 = 'integrated_snn_res.0.8',
                   col2 = 'annotationsingleRConsensusCell')
give_table_matches(combined_dataset_light_df = combined_dataset_light[[]],
                   col1 = 'integrated_snn_res.0.8', col2 = 'annotationsingleRConsensusCell',
                   prune=T)
give_table_matches(combined_dataset_light_df = combined_dataset_lightCharly[[]], col1 = 'integrated_snn_res.4',
                   col2 = 'annotationsingleRConsensusCell', prune = T, prune_val=0.3)


clustree(combined_dataset_light)
# clustree(combined_dataset_lightCharly, prefix = )

combined_dataset_light_clustree <- combined_dataset_light
combined_dataset_light_clustree$clustreecomparison1 <- combined_dataset_light_clustree$integrated_snn_res.0.8
combined_dataset_light_clustree$clustreecomparison2 <- combined_dataset_light_clustree$annotationsingleRConsensusCell
## remove cells with NAs
combined_dataset_light_clustree <- subset(combined_dataset_light_clustree, cells=dimnames(combined_dataset_light_clustree)[[2]][-which(is.na(combined_dataset_light_clustree$clustreecomparison1) | is.na(combined_dataset_light_clustree$clustreecomparison2))])

combined_dataset_light_clustree ## 81837 samples
combined_dataset_light ## 82584 samples

## I could remove a lot of data from <combined_dataset_light_clustree> to make it a lighter object

clustree(combined_dataset_light_clustree, prefix = 'clustreecomparison')

riverplot(list(nodes=cbind.data.frame(x=combined_dataset_light_clustree$integrated_snn_res.0.8,
                           y=combined_dataset_light_clustree$annotationsingleRConsensusCell)))

## this hasn't worked
## https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
ggplot(plyr::ddply(cbind.data.frame(x=combined_dataset_light_clustree$integrated_snn_res.0.8,
                                    y=combined_dataset_light_clustree$annotationsingleRConsensusCell),.(x,y),nrow),
       aes(y = V1, axis1 = x, axis2 = y)) +
  geom_alluvium(aes(fill = x), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey")+
  geom_label(stat = "stratum", aes(label = after_stat(stratum)))

