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

##-----------------------------------------------------------------------------------------------------------------##
rename_observations <- function(i){
  splt <- strsplit(i, '-')[[1]]
  paste0(gsub(".RDS", "", gsub("clusters_nPCs", "", basename(input_fles_clusters)))[as.numeric(splt[1])],
         '-',
         gsub("integrated_snn_res.", "", colnames(files_clusters[[1]]))[as.numeric(splt[2])])
}

reorder_with_colnames <- function(m, reorder){
  m <- m[,reorder]
  colnames(m) <- colnames(m)[reorder]
  m
}

add_symmetry_and_dcast <- function(m){
  m <- m[!is.na(m$X3),]
  m <- rbind(m, reorder_with_colnames(m, c(2,1,3)))
  m <- m %>% distinct %>% data.frame
  ## make symmetrical
  m_dcast <- dcast(m, X1~X2, value.var = "X3")
  rownames(m_dcast) <- m_dcast[,1]
  m_dcast[,1] <- NULL
  return(m_dcast)
}

colname_as_numeric <- function(i){
  rw_i <- rownames(i)
  x <- apply(i, 2, as.numeric)
  rownames(x) <- rw_i
  x
}

give_complexheatmap_with_anno <- function(m, title_arg='value', both_anno=F, no_anno=F, ...){
  require(ComplexHeatmap)
  m_meta <- list(rowmeta=data.frame(row.names = rownames(m),
                                    nPCs = as.numeric(sapply(rownames(m), function(i) strsplit(i, '-')[[1]][1])),
                                    resolution = as.numeric(sapply(rownames(m), function(i) strsplit(i, '-')[[1]][2])) ),
                 colmeta=data.frame(row.names = colnames(m),
                                    nPCs = as.numeric(sapply(colnames(m), function(i) strsplit(i, '-')[[1]][1])),
                                    resolution = as.numeric(sapply(colnames(m), function(i) strsplit(i, '-')[[1]][2])) ))
  
  if(no_anno){
    ComplexHeatmap::Heatmap(m, heatmap_legend_param = list(title = title_arg),
                            ...)
  }else{
    if(both_anno){
      ComplexHeatmap::Heatmap(m, heatmap_legend_param = list(title = title_arg),
                              top_annotation = HeatmapAnnotation(df = m_meta$colmeta),
                              left_annotation = rowAnnotation(df = m_meta$rowmeta),
                              ...)
    }else{
      ComplexHeatmap::Heatmap(m, heatmap_legend_param = list(title = title_arg),
                              top_annotation = HeatmapAnnotation(df = m_meta$colmeta),
                              ...)
    }
  }
}

##-----------------------------------------------------------------------------------------------------------------##

##-----------------------------------------------------------------------------------------------------------------##

##' tile of nPCs (x) and resolution (y) showing the similarity in clustering between two methods. this could be assessed by
##' several metrics of clustering
##' https://stats.stackexchange.com/questions/95782/what-are-the-most-common-metrics-for-comparing-two-clustering-algorithms-especi


#calculate Rand index between clustering methods
subset_cells <- sample(1:nrow(files_clusters[[1]]), size = 2000)

nPC_res_combination <- expand.grid(1:length(files_clusters), 1:ncol(files_clusters[[1]]))
nPC_res_combination_paste <- apply(nPC_res_combination, 1, paste0, collapse= '-')
# nPC_res_combination_paste <- factor(nPC_res_combination_paste)
randindices_subset <- data.frame(matrix(NA, nrow=length(nPC_res_combination_paste)*length(nPC_res_combination_paste), ncol=3))
count <- 1
for(iidx in 1:length(nPC_res_combination_paste)){
  for(jidx in 1:length(nPC_res_combination_paste)){
    cat(i, '\t', j, '\n')
    if((iidx)>(jidx)){ ## so that it's non-redundant
      i <- nPC_res_combination_paste[iidx]
      j <- nPC_res_combination_paste[jidx]
      isplit <- as.numeric(strsplit(as.character(i), '-')[[1]])
      jsplit <- as.numeric(strsplit(as.character(j), '-')[[1]])
      .x <- rand.index(as.numeric(files_clusters[[isplit[1]]][subset_cells,isplit[2]]),
                       as.numeric(files_clusters[[jsplit[1]]][subset_cells,jsplit[2]]))
    }else{
      .x <- NA
    }
    randindices_subset[count,] <- c(i,j,.x)
    count <- count+1
  }
}
randindices_subset[,1] <- sapply(randindices_subset[,1], rename_observations)
randindices_subset[,2] <- sapply(randindices_subset[,2], rename_observations)

# saveRDS(randindices_subset, paste0(robject_folder, 'randindices_subset.RDS'))
randindices_subset <- readRDS(paste0(robject_folder, 'randindices_subset.RDS'))


randindices_subset_dcast <- add_symmetry_and_dcast(randindices_subset)
randindices_subset_dcast <- colname_as_numeric(randindices_subset_dcast)


image(as(randindices_subset_dcast, 'matrix'))

# randindices_subset_dcast <- randindices_subset_dcast[order(rowSums(is.na(randindices_subset_dcast)), decreasing = T),]
# randindices_subset_dcast <- randindices_subset_dcast[,order(colSums(is.na(randindices_subset_dcast)))]
## make symmetrical
# for(i in 1:nrow(randindices_subset_dcast)){
#   for(j in i:ncol(randindices_subset_dcast)){
#     randindices_subset_dcast[i,j] <- randindices_subset_dcast[j,i]
#   }
# }

pheatmap::pheatmap(randindices_subset_dcast)
pheatmap::pheatmap(randindices_subset_dcast, cluster_rows = F, cluster_cols = F)
randindices_subset_dcast

# res_vect <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5)

randindices_subset_dcast_meta <- list(rowmeta=data.frame(row.names = rownames(randindices_subset_dcast),
                                                         nPCs = as.numeric(sapply(rownames(randindices_subset_dcast), function(i) strsplit(i, '-')[[1]][1])),
                                                         resolution = (sapply(rownames(randindices_subset_dcast), function(i) strsplit(i, '-')[[1]][2])) ),
                                      colmeta=data.frame(row.names = colnames(randindices_subset_dcast),
                                                         nPCs = as.numeric(sapply(colnames(randindices_subset_dcast), function(i) strsplit(i, '-')[[1]][1])),
                                                         resolution = as.numeric(sapply(colnames(randindices_subset_dcast), function(i) strsplit(i, '-')[[1]][2])) ))

image(as(randindices_subset_dcast, 'matrix'))
## add Aerts and any other independent classification
try(dev.off())
pdf(paste0(folder_results, input_objs, "/heatmap_clusters_PCs_resolution_randindices.pdf"), height = 20, width = 16)
ComplexHeatmap::Heatmap(randindices_subset_dcast,
                        top_annotation = HeatmapAnnotation(df = randindices_subset_dcast_meta$colmeta),
                        left_annotation = rowAnnotation(df = randindices_subset_dcast_meta$rowmeta))
dev.off()

ComplexHeatmap::Heatmap(randindices_subset_dcast,
                        top_annotation = HeatmapAnnotation(df = randindices_subset_dcast_meta$colmeta),
                        left_annotation = rowAnnotation(df = randindices_subset_dcast_meta$rowmeta),
                        cluster_rows = F, cluster_columns = F) ## good clustering

image(table(files_clusters[[1]][,1:2]))

rand.index(files_clusters[[1]][1:2000,1], files_clusters[[1]][1:2000,2])

matrix_subset_grep <- function(m, string){
  m[grepl(string, rownames(m)),grepl(string, colnames(m))]
}

give_complexheatmap_with_anno(matrix_subset_grep(randindices_subset_dcast, '^10-'))

##---------------


VIdists <- data.frame(matrix(NA, nrow=length(nPC_res_combination_paste)*length(nPC_res_combination_paste), ncol=3))
count <- 1
for(iidx in 1:length(nPC_res_combination_paste)){
  for(jidx in 1:length(nPC_res_combination_paste)){
    cat(i, '\t', j, '\n')
    if((iidx)>(jidx)){ ## so that it's non-redundant
      i <- nPC_res_combination_paste[iidx]
      j <- nPC_res_combination_paste[jidx]
      isplit <- as.numeric(strsplit(as.character(i), '-')[[1]])
      jsplit <- as.numeric(strsplit(as.character(j), '-')[[1]])
      .x <- mcclust::vi.dist(as.numeric(files_clusters[[isplit[1]]][,isplit[2]]),
                             as.numeric(files_clusters[[jsplit[1]]][,jsplit[2]]))
    }else{
      .x <- NA
    }
    VIdists[count,] <- c(i,j,.x)
    count <- count+1
  }
}

# saveRDS(VIdists, paste0(robject_folder, 'VIdists.RDS'))
VIdists <- readRDS(paste0(robject_folder, 'VIdists.RDS'))
VIdists[,1] <- sapply(VIdists[,1], rename_observations)
VIdists[,2] <- sapply(VIdists[,2], rename_observations)

VIdists_dcast <- add_symmetry_and_dcast(VIdists)
VIdists_dcast <- colname_as_numeric(VIdists_dcast)

give_complexheatmap_with_anno(matrix_subset_grep(VIdists_dcast, '^39-'))
give_complexheatmap_with_anno(matrix_subset_grep(VIdists_dcast, '-9$'))
give_complexheatmap_with_anno(VIdists_dcast)
give_complexheatmap_with_anno(VIdists_dcast, cluster_rows = F, cluster_columns = F)
colnames(VIdists_dcast)

stopifnot(all(colnames(VIdists_dcast) == colnames(randindices_subset_dcast)))
stopifnot(all(rownames(VIdists_dcast) == rownames(randindices_subset_dcast)))

plot(unlist(VIdists_dcast), unlist(randindices_subset_dcast))

ggplot(data.frame(VIdists_dcast=as.vector(unlist(VIdists_dcast)),
                  randindices_subset_dcast=as.vector(unlist(randindices_subset_dcast))),
       aes(VIdists_dcast, randindices_subset_dcast)) +
  # geom_hex(bins = 100) +
  geom_pointdensity() +
  scale_color_scico(palette = "devon", direction = -1, end = 0.9)+
  theme_bw()+labs(x='VI distance', y='Rand index (subset of cells)')
ggsave(paste0(folder_results, input_objs, "/clusters_PCs_resolution_Randindex_VIdist_scatter.pdf"), height = 5, width = 6)

par(mfrow=c(1,2))
image(as(randindices_subset_dcast, 'matrix'))
image(as(VIdists_dcast, 'matrix'))

sort_by_res <- function(i){
  i <- i[,gtools::mixedorder(sapply(colnames(i), function(i) paste0(rev(strsplit(i, '-')[[1]]), collapse='-')))]
  i <- i[gtools::mixedorder(sapply(rownames(i), function(i) paste0(rev(strsplit(i, '-')[[1]]), collapse='-'))),]
  i
}

pdf(paste0(folder_results, input_objs, "/heatmap_clusters_PCs_resolution_sorted_VIdist.pdf"), height = 20, width = 16)
give_complexheatmap_with_anno(sort_by_res(VIdists_dcast), cluster_rows = F, cluster_columns = F)
dev.off()

pdf(paste0(folder_results, input_objs, "/heatmap_clusters_PCs_resolution_sorted_randindices.pdf"), height = 20, width = 16)
give_complexheatmap_with_anno(sort_by_res(randindices_subset_dcast), cluster_rows = F, cluster_columns = F)
dev.off()

pdf(paste0(folder_results, input_objs, "/heatmap_clusters_PCs_resolution_VIdist.pdf"), height = 20, width = 16)
give_complexheatmap_with_anno(VIdists_dcast, cluster_rows = T, cluster_columns = T)
dev.off()

pdf(paste0(folder_results, input_objs, "/heatmap_clusters_PCs_resolution_randindices.pdf"), height = 20, width = 16)
give_complexheatmap_with_anno(randindices_subset_dcast, cluster_rows = T, cluster_columns = T)
dev.off()

### get the number of clusters in each of the combinations

num_clusters <- outer(1:length(files_clusters), 1:length(files_clusters[[1]]),
                      Vectorize(function(i,j){
                        length(unique(files_clusters[[i]][[j]]))
                      }
                      ))
rownames(num_clusters) <- gsub(".RDS", "", gsub("clusters_nPCs", "", basename(input_fles_clusters)))
colnames(num_clusters) <- gsub("integrated_snn_res.", "", colnames(files_clusters[[1]]))
num_clusters
num_clusters <- num_clusters[,order(colnames(num_clusters))]
num_clusters <- num_clusters[order(rownames(num_clusters)),]

pdf(paste0(folder_results, input_objs, "/heatmap_clusters_PCs_resolution_sorted_numclusters.pdf"), height = 4, width = 5)
give_complexheatmap_with_anno(num_clusters, cluster_rows = F, cluster_columns = F, title='Number of clusters', no_anno=T)
dev.off()

pheatmap(num_clusters)
