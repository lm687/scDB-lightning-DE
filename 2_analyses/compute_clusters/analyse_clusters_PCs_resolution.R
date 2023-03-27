library(fossil)
library(mcclust)
robject_folder <- "~/projects/lighting/data/robjects/" ## CCB cluster
genome <- 'dmel649ChrimsonV2'
input_objs <- paste0('KB-', genome)
robject_folder <- paste0(robject_folder, input_objs, '/')
folder_results <- "~/projects/lighting/3_results/" ## CCB cluster

input_fles_clusters <- list.files("/home/l/lmorrill/projects/lighting/data/robjects/KB-dmel649ChrimsonV2/clustering_nPCs/", full.names = T)

files_clusters <- lapply(input_fles_clusters, readRDS)


##' tile of nPCs (x) and resolution (y) showing the similarity in clustering between two methods. this could be assessed by
##' several metrics of clustering
##' https://stats.stackexchange.com/questions/95782/what-are-the-most-common-metrics-for-comparing-two-clustering-algorithms-especi


#calculate Rand index between clustering methods
subset_cells <- sample(1:nrow(files_clusters[[1]]), size = 2000)

nPC_res_combination <- expand.grid(1:length(files_clusters), 1:ncol(files_clusters[[1]]))
nPC_res_combination_paste <- apply(nPC_res_combination, 1, paste0, collapse= '-')
nPC_res_combination_paste <- factor(nPC_res_combination_paste)
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

randindices_subset[,1] <- as.character(nPC_res_combination_paste[randindices_subset[,1]])
randindices_subset[,2] <- as.character(nPC_res_combination_paste[randindices_subset[,2]])
# saveRDS(randindices_subset, paste0(robject_folder, 'randindices_subset.RDS'))
randindices_subset <- readRDS(paste0(robject_folder, 'randindices_subset.RDS'))


reorder_with_colnames <- function(m, reorder){
  m <- m[,reorder]
  colnames(m) <- colnames(m)[reorder]
  m
}

add_symmetry_and_dcast <- function(m){
  m <- randindices_subset[!is.na(randindices_subset$X3),]
  m <- rbind(m, reorder_with_colnames(m, c(2,1,3)))
  m <- m %>% distinct %>% data.frame
  ## make symmetrical
  m_dcast <- dcast(m, X1~X2, value.var = "X3")
  rownames(m_dcast) <- m_dcast[,1]
  m_dcast[,1] <- NULL
  return(m_dcast)
}

randindices_subset_dcast <- add_symmetry_and_dcast(randindices_subset_noNA)


randindices_subset_dcast['9-9', '1-2']
randindices_subset_dcast['1-2', '9-9']

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

res_vect <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5)

randindices_subset_dcast_meta <- list(rowmeta=data.frame(row.names = rownames(randindices_subset_dcast),
                                                         nPCs = as.numeric(sapply(rownames(randindices_subset_dcast), function(i) strsplit(i, '-')[[1]][1])),
                                                         resolution = res_vect[as.numeric(sapply(rownames(randindices_subset_dcast), function(i) strsplit(i, '-')[[1]][2]))] ),
                                      colmeta=data.frame(row.names = colnames(randindices_subset_dcast),
                                                         nPCs = as.numeric(sapply(colnames(randindices_subset_dcast), function(i) strsplit(i, '-')[[1]][1])),
                                                         resolution = res_vect[as.numeric(sapply(colnames(randindices_subset_dcast), function(i) strsplit(i, '-')[[1]][2]))] ))

library(ComplexHeatmap)
## add Aerts and any other independent classification
try(dev.off())
pdf(paste0(folder_results, input_objs, "/heatmap_clusters_PCs_resolution.pdf"), height = 20, width = 16)
ComplexHeatmap::Heatmap(randindices_subset_dcast,
                        top_annotation = HeatmapAnnotation(df = randindices_subset_dcast_meta$colmeta),
                        left_annotation = rowAnnotation(df = randindices_subset_dcast_meta$rowmeta))
dev.off()

give_complexheatmap_with_anno <- function(m){
  require(ComplexHeatmap)
  m_meta <- list(rowmeta=data.frame(row.names = rownames(m),
                                                           nPCs = as.numeric(sapply(rownames(m), function(i) strsplit(i, '-')[[1]][1])),
                                                           resolution = res_vect[as.numeric(sapply(rownames(m), function(i) strsplit(i, '-')[[1]][2]))] ),
                                        colmeta=data.frame(row.names = colnames(m),
                                                           nPCs = as.numeric(sapply(colnames(m), function(i) strsplit(i, '-')[[1]][1])),
                                                           resolution = res_vect[as.numeric(sapply(colnames(m), function(i) strsplit(i, '-')[[1]][2]))] ))
  
  ComplexHeatmap::Heatmap(m,
                          # top_annotation = HeatmapAnnotation(df = m_meta$colmeta),
                          left_annotation = rowAnnotation(df = m_meta$rowmeta))
}

##I DON'T THINK THIS CAN BE CORRECT - PROBABLY THE RESOLUTIONS ARE NOT CORRECT ?????


image(table(files_clusters[[1]][,1:2]))

rand.index(files_clusters[[1]][1:2000,1], files_clusters[[1]][1:2000,2])

matrix_subset_grep <- function(m, string){
  m[grepl(string, rownames(m)),grepl(string, colnames(m))]
}

give_complexheatmap_with_anno(matrix_subset_grep(randindices_subset_dcast, '^2-'))

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

VIdists[,1] <- as.character(nPC_res_combination_paste[VIdists[,1]])
VIdists[,2] <- as.character(nPC_res_combination_paste[VIdists[,2]])
# saveRDS(VIdists, paste0(robject_folder, 'VIdists.RDS'))
VIdists_dcast <- add_symmetry_and_dcast(VIdists)

