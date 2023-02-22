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
  # seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:60, umap.method='umap-learn') ## couldn't install it locally (??)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:60, umap.method='uwot')
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

give_top_logf_genes_val <- function(i, include_separation_row=T, n=6L){
  cat('Sorting by p-value\n')
  i <- i[order(i$p_val_adj, na.last = T),]
  return(head(i, n=n))
}


give_cluster_specific_DE <- function(seurat_obj, cluster_name=NULL, group.by_arg=NULL, ident.1_arg, ident.2_arg){
  warning('Function <FoldChanges> has been replaced by <FindMarkers> in order to find p-values and test for DE')
  stopifnot(!is.null(cluster_name))
  stopifnot(!is.null(group.by_arg))
  clusters <- sort(unique(seurat_obj@meta.data[,cluster_name]))
  res <- lapply(clusters, function(cluster_it){
    subset_cells <- subset(seurat_obj, cells = which(seurat_obj@meta.data[,cluster_name] == cluster_it))
    Seurat::FindMarkers(subset_cells, ident.1=ident.1_arg, ident.2_arg, group.by=group.by_arg)
  })
  names(res) <- clusters
  return(res)
}

normalise_rw <- function(i){
  if(is.null(dim(i))){
    i/sum(i)
  }else{
    sweep(i, 1, rowSums(i), '/')
  }
}


normalise_cl <- function(i){
  sweep(i, 2, colSums(i), '/')
}

rownames_to_col <- function(i){
  data.frame(names=rownames(i), i)
}

firstcol_to_rownames <- function(i){
  rownames(i) <- i[,1]
  i[,-1]
}

na_to_zeros <- function(i){
  i[is.na(i)] <- 0
  i
}

relevel_by_value_column <- function(i){
  i$names <- factor(i$names, levels=i$names[order(i$value)])
  return(i)
}

give_name_conversion_file <- function(return_from_fb=T,
                                      path_tsv_genes="/home/l/lmorrill/projects/general/fb_synonym_fb_2022_06_cut.tsv"){
  # library("AnnotationDbi")
  # library("org.Dm.eg.db")
  
  
  if(return_from_fb){
    # fly_genes <- read.table("/home/l/lmorrill/projects/general/fb_synonym_fb_2022_06_cut.tsv", sep = "\t", comment.char = "^#")
    # fly_genes <- read.table("/home/l/lmorrill/projects/general/fb_synonym_fb_2022_06_cut.tsv", sep = "\t")
    fly_genes <- read.table(path_tsv_genes)
    # fly_genes[grepl('fly_genes', fly_genes$V1)]
    # dim(fly_genes)
    # system("wc -l /home/l/lmorrill/projects/general/fb_synonym_fb_2022_06.tsv")
    ## rows are missing
    
    # fly_genes[1:2,1]
    # fly_genes[,2][match('FBgn0024733', fly_genes[,1])]
    # fly_genes[,2][match('FBgn0031081', fly_genes[,1])]

    return(fly_genes)
   
    ## there are many NAs in gene names
  }else{
    
    # https://www.biostars.org/p/70821/
    # resadj <- as.data.frame(resadj)
    # res$symbol <- mapIds(org.Dm.eg.db, 
    #                      keys=row.names(res), 
    #                      column="SYMBOL", 
    #                      keytype="ENSEMBL",
    #                      multiVals="first")
    # 
    # res$entrez <- mapIds(org.Dm.eg.db, 
    #                      keys=row.names(res), 
    #                      column="ENTREZID", 
    #                      keytype="ENSEMBL",
    #                      multiVals="first")
    # 
    # res$name =   mapIds(org.Dm.eg.db,
    #                     keys=row.names(res), 
    #                     column="GENENAME",
    #                     keytype="ENSEMBL",
    #                     multiVals="first")
    # write.csv(res, file = "/home/l/lmorrill/projects/general/results_FlyBaseIDS_GeneNames.csv", sep = "\t", quote = F)
  }
}

convert_FB_to_name <- function(i){
  name_conversion_file$V2[match(i, name_conversion_file$V1)]
}


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
  a+theme_bw()
}

add_genenames <- function(i){
  lapply(i, function(j) cbind.data.frame(gene=rownames(j), j))
}