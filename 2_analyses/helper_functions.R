## Functions
splitUMAPPlot <- function(seurat_object, group.by=NULL, dim_red='umap', colour_by_feature=NULL){
  stopifnot(!is.null(group.by))
  dimred_mat <- cbind.data.frame(seurat_object@reductions[[dim_red]]@cell.embeddings,
                                 group=seurat_object[[group.by]])
  if(!is.null(colour_by_feature)){
    dimred_mat$colfill <- as.vector(seurat_object@assays$RNA[colour_by_feature])
    colnames(dimred_mat) <- c('dim1', 'dim2', 'group', 'colfill')
    
  }else{
    colnames(dimred_mat) <- c('dim1', 'dim2', 'group')
  }
  dimred_mat$col=T
  
  extra_background <- do.call('rbind', replicate(length(unique(dimred_mat$group)), dimred_mat, simplify = F))
  extra_background$col <- F
  extra_background$group2 <- extra_background$group
  extra_background$group <- rep(unique(dimred_mat$group), each=nrow(dimred_mat))
  extra_background <- extra_background[!(extra_background$group == extra_background$group2),]
  extra_background$group2 <- NULL
  dimred_mat <- rbind(dimred_mat, extra_background)
  dimred_mat <- dimred_mat[order(dimred_mat$group),]
  # dimred_mat$col <- factor(dimred_mat$col, levels=c('TRUE', 'FALSE'))
  
  if(is.null(colour_by_feature)){
    ggplot(dimred_mat[dimred_mat$col,],
           aes(x=dim1, y=dim2, col=col))+
      geom_point(size=0.3, data=dimred_mat[!dimred_mat$col,], alpha=0.2)+
      geom_point(shape=1)+
      facet_wrap(.~group)+theme_bw()+labs(col='GE')
  }else{
    colour_by_feature
    ggplot(dimred_mat[dimred_mat$col,],
           aes(x=dim1, y=dim2, col=colfill))+
      geom_point(size=0.3, data=dimred_mat[!dimred_mat$col,], alpha=0.2)+
      geom_point(shape=1)+
      facet_wrap(.~group)+theme_bw()+labs(col='GE')
  }
}

add_dimred <- function(seurat_obj){
  cat('Scaling...\n')
  seurat_obj <- ScaleData(seurat_obj)
  cat('Running PCA...\n')
  seurat_obj <- RunPCA(seurat_obj, npcs = 60)
  cat('Running UMAP\n')
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:60)
  cat('Finding neighbours...\n')
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:60)
  cat('Finding clusters...\n')
  seurat_obj <- FindClusters(seurat_obj, resolution = 4)
  cat('Running t-SNE...\n')
  seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = 1:60)
  return(seurat_obj)
}


give_top_logf_genes <- function(i, include_separation_row=T, n=6L){
  i <- i[order(i[,1]),]
  if(include_separation_row){
    return(rbind(head(i, n=n), rep('...', ncol(i)), tail(i, n=n)))
  }else{
    return(rbind(head(i, n=n), tail(i, n=n)))
  }
}

give_cluster_specific_DE <- function(seurat_obj, cluster_name=NULL, group.by_arg=NULL, ident.1_arg, ident.2_arg){
  stopifnot(!is.null(cluster_name))
  stopifnot(!is.null(group.by_arg))
  clusters <- sort(unique(Chrimsonposcells@meta.data[,cluster_name]))
  res <- lapply(clusters, function(cluster_it){
    subset_cells <- subset(Chrimsonposcells, cells = which(Chrimsonposcells@meta.data[,cluster_name] == cluster_it))
    Seurat::FoldChange(subset_cells, ident.1=ident.1_arg, ident.2_arg, group.by=group.by_arg)
  })
  names(res) <- clusters
  return(res)
}