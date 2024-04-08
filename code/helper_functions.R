## Functions
splitUMAPPlot <- function(...){
  warning('Use splitDimRedPlot instead ')
}
splitDimRedPlot <- function(seurat_object, group.by=NULL, dim_red='umap', colour_by_feature=NULL, plot_only_true=F, show_other=T, plot_legend=T){
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
  
  if(!show_other){
    ## don't plot the first group (which should be the background) as a facet
    if(any(dimred_mat$group == 'Other')){
      dimred_mat <- dimred_mat[dimred_mat$group != 'Other',]
    }else{
      dimred_mat <- dimred_mat[dimred_mat$group == unique(dimred_mat$group)[1],]
    }
  }
  
  if(plot_only_true){
    dimred_mat <- dimred_mat[dimred_mat$col,]
  }
  if(is.null(colour_by_feature)){
    resplot <- ggplot(dimred_mat[dimred_mat$col,],
           aes(x=dim1, y=dim2, col=col))+
      geom_point(size=0.3, data=dimred_mat[!dimred_mat$col,], alpha=0.2)+
      # geom_point(shape=0)+
      geom_point()+
      facet_wrap(.~group)+theme_bw()+labs(col='GE')+scale_color_manual(values = c('#e8e8e8', '#00cd9e')) ## grey, colour
  }else{
    colour_by_feature
    resplot <- ggplot(dimred_mat[dimred_mat$col,],
           aes(x=dim1, y=dim2, col=colfill))+
      geom_point(size=0.3, data=dimred_mat[!dimred_mat$col,], alpha=0.2)+
      geom_point(shape=1)+
      facet_wrap(.~group)+theme_bw()+labs(col='GE')
  }
  if(!plot_legend){
    resplot <- resplot+guides(col='none')
  }
  resplot
}

add_dimred <- function(seurat_obj, nPCs=60, resolution_clusters=4, compute_neighbours=T, scale_and_PCA=T){
  if(scale_and_PCA){
    cat('Scaling...\n')
    seurat_obj <- ScaleData(seurat_obj)
    cat('Running PCA...\n')
    seurat_obj <- RunPCA(seurat_obj, npcs = 60)
  }
  cat('Running UMAP\n')
  # seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:60, umap.method='umap-learn') ## couldn't install it locally (??)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = 1:nPCs, umap.method='uwot')
  if(compute_neighbours){
    cat('Finding neighbours...\n')
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:nPCs)
    cat('Finding clusters...\n')
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution_clusters)
    cat('Running t-SNE...\n')
    seurat_obj <- RunTSNE(seurat_obj, reduction = "pca", dims = 1:nPCs)
  }
  return(seurat_obj)
}


rename_metadata_column <- function(seurat_obj, previous_colname, current_colname){
  if(previous_colname %in% colnames(seurat_obj@meta.data)){
    seurat_obj@meta.data[,current_colname] <- seurat_obj@meta.data[,previous_colname]
    seurat_obj@meta.data[,previous_colname] <- NULL
  }else{
    warning('Column not found!')
  }
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
  # cat('Here\n')
  res <- lapply(clusters, function(cluster_it){
    idx_select <- which(seurat_obj@meta.data[,cluster_name] == cluster_it)
    # print(idx_select)
    subset_cells <- subset(seurat_obj, cells = idx_select)
    cat('Number of cells in each DE group: ', table(subset_cells@meta.data[,group.by_arg]), '\n')
    if( (length(unique(subset_cells@meta.data[,group.by_arg])) == 1) | any(table(subset_cells$stim) < 3) ){
      
      if(length(unique(subset_cells@meta.data[,group.by_arg])) == 1)     warning('There is perfect separation between DE groups in cluster', cluster_it, '\n')
      if(any(table(subset_cells$stim) < 3))    warning('There are not enough cells in each condition group in cluster', cluster_it, '\n')
      
      ## there is perfect separation or another problem - return empty dataframe
      setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"))
    }else{
      Seurat::FindMarkers(subset_cells, ident.1=ident.1_arg, ident.2_arg, group.by=group.by_arg)
    }
  })
  names(res) <- clusters
  cat('Returning cluster-specific DEGs\n')
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

remove_cols_with_na <- function(i){
  i[,!is.na(colSums(i))]
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
                                      path_tsv_genes="/project/sims-lab/lmorrill/data/general/genes/fb_synonym_fb_2022_06_cut.tsv"){
  # library("AnnotationDbi")
  # library("org.Dm.eg.db")
  
  cat('The path to the transcript to gene conversion file is ', path_tsv_genes, '\n')
  
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

rownames_as_col <- function(i){
  cbind.data.frame(name=rownames(i), i)
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

remove_na <- function(i) i[!is.na(i)]

umap_single_facet_with_topmarker <- function(seurat_obj, markers, seurat_name_clusters, density_map=F){
  require(ggrepel)
  reduction_df <- data.frame(Reductions(seurat_obj, 'umap')@cell.embeddings)
  reduction_df$clusters <- seurat_obj[[seurat_name_clusters]][,1]
  reduction_df$topmarker <- NA
  if(!is.null(levels(markers$cluster))){
    ## it's a factor
    markers$cluster <- as.character(markers$cluster)
  }
  if(!is.null(levels(reduction_df$clusters))){
    ## it's a factor
    reduction_df$clusters <- as.character(reduction_df$clusters)
  }
  
  reduction_df <- rbind(reduction_df,
                        cbind.data.frame(UMAP_1 = sapply(unique(reduction_df$clusters), function(i) mean(reduction_df$UMAP_1[which(reduction_df$cluster == i)])),
                                         UMAP_2= sapply(unique(reduction_df$clusters), function(i) mean(reduction_df$UMAP_2[which(reduction_df$cluster == i)])),
                                         clusters = unique(reduction_df$clusters),
                                         topmarker = markers$gene[sapply(unique(reduction_df$clusters), function(i) which(markers$cluster == i)[1])])
  )
  
  warning('Non latin ASCII characters have been changed')
  reduction_df$topmarker <- stringi::stri_trans_general(reduction_df$topmarker, "latin-ascii")
  
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


give_DE_analysis_lightON_lightOFF <- function(dataset, dataset_name, name_clusters=NULL, add_name='', FoldChange=T){
  require(ggridges)
  
  bl <- colorRampPalette(c("navy","royalblue","lightskyblue"))(200)                      
  re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
  
  if(is.null(dataset@meta.data$stim)){
    stop('dataset should contain column <stim> in metadata')
  }
  dataset@meta.data$stim_light <- (dataset@meta.data$stim %in% c('G1', 'G3'))
  print(table(dataset@meta.data$stim_light))
  dataset@meta.data$stim_light[dataset@meta.data$stim_light] <- 'Lights on'
  dataset@meta.data$stim_light[dataset@meta.data$stim_light == 'FALSE'] <- 'Lights off'
  if(!is.null(name_clusters)) dataset$seurat_clusters <- dataset[[name_clusters]][,1]
  if(FoldChange){
    # FoldChange: on presence/absence
    DE_Chrimsonposcells_lights <- Seurat::FoldChange(dataset, ident.1='Lights on', ident.2='Lights off', group.by='stim_light')
  }else{
    # FindMarkers: on expression
    DE_Chrimsonposcells_lights <- Seurat::FindMarkers(dataset, ident.1='Lights on', ident.2='Lights off', group.by='stim_light')
  }
  DE_Chrimsonposcells_lights
  DE_Chrimsonposcells_lights_topgenes <- give_top_logf_genes(DE_Chrimsonposcells_lights)
  xtable::xtable(DE_Chrimsonposcells_lights_topgenes)
  
  system(paste0("mkdir -p ", folder_results, input_objs, "/", dataset_name, add_name))
  
  cat('Saving ', paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_featureplot.png"), '\n')
  FeaturePlot(dataset, features = rownames(DE_Chrimsonposcells_lights_topgenes))
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_featureplot.png"), height = 6, width = 8)
  
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
  if(local){
    saveRDS(DE_Chrimsonposcells_lights_per_cluster, paste0("/Users/lenamorrill/Documents/projects/lightning/github-repo-lightning-DE/3_results_local/objects/", input_objs, "/", dataset_name, add_name, "_DEgenes_per_cluster.RDS"))
    # DE_Chrimsonposcells_lights_per_cluster <- readRDS(paste0("/Users/lenamorrill/Documents/projects/lightning/github-repo-lightning-DE/3_results_local/objects/", input_objs, "/", dataset_name, "_DEgenes_per_cluster.RDS"))
  }else{
    cat('Saving ', paste0("~/projects/lighting/data/robjects/", input_objs, "/", dataset_name, add_name, "_DEgenes_per_cluster.RDS"), '\n')
    saveRDS(DE_Chrimsonposcells_lights_per_cluster, paste0("~/projects/lighting/data/robjects/", input_objs, "/", dataset_name, add_name, "_DEgenes_per_cluster.RDS"))
  }
  
  cat('Selecting clusters without perfect separation in the tested condition\n')
  cat(sum(sapply(DE_Chrimsonposcells_lights_per_cluster, nrow) == 0), 'cluster/s had perfect separation and was/were removed.\n')
  DE_Chrimsonposcells_lights_per_cluster <- DE_Chrimsonposcells_lights_per_cluster[sapply(DE_Chrimsonposcells_lights_per_cluster, nrow) > 0]
  
  cat('Selecting DEGs\n')
  DE_Chrimsonposcells_lights_per_cluster_topgenes <- lapply(DE_Chrimsonposcells_lights_per_cluster,function(i) i[i$p_val_adj <= 0.05,])
  cat('Selexting top DEGs\n')
  DE_Chrimsonposcells_lights_per_cluster_topgenes <- lapply(DE_Chrimsonposcells_lights_per_cluster_topgenes, function(i) data.frame(gene=rownames(i), avg_log2FC=i[,'avg_log2FC']))
  DE_Chrimsonposcells_lights_per_cluster_topgenes <- melt(DE_Chrimsonposcells_lights_per_cluster_topgenes)
  table(DE_Chrimsonposcells_lights_per_cluster_topgenes$variable)
  # print('here')
  cat('Dcasting DEGs\n')
  DE_Chrimsonposcells_lights_per_cluster_topgenesdcast <- dcast(DE_Chrimsonposcells_lights_per_cluster_topgenes, L1~gene, value.var = "value")
  # print('here 2')
  DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[is.na(DE_Chrimsonposcells_lights_per_cluster_topgenesdcast)] <- 0 ## semi-controversial
  # print('here 3')
  
  DE_Chrimsonposcells_lights_per_cluster_topgenes$gene <- factor(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene, levels=names(sort(table(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene))))
  ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes, aes(y=gene, x=L1, fill=value))+
    scale_fill_gradientn(colours=c(bl,"white", re), na.value = "grey98")+
    geom_tile()+ggtitle('Differentially expressed genes between lights ON/OFF per cell cluster')+
    theme_bw()+labs(x='Cluster', y='Gene')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_percluster_DEgenes.png"),
         height = min(length(levels(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene)), 15), width = 8)
  
  
  ## to find clustering
  select_two_or_more_active <- function(i){
    i[,colSums(i>0) > 2]
  }
  
  try({
    try(dev.off())
    pdf(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_percluster_DEgenes_subsetgenes_cluster.pdf"), height = 3, width = 4)
    # print(pheatmap::pheatmap(t(as(select_two_or_more_active(DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1]), 'matrix'))))
    print(pheatmap::pheatmap(t(as((DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1]), 'matrix'))))
    dev.off()
    # print('here 4')
  })
  
  try({
    try(dev.off())
    pdf(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_percluster_DEgenes_subsetgenes_cluster_large.pdf"),
        height = min(length(levels(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene)), 15), width = 8)
    # print(pheatmap::pheatmap(t(as(select_two_or_more_active(DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1]), 'matrix'))))
    print(pheatmap::pheatmap(t(as((DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1]), 'matrix'))))
    dev.off()
    # print('here 4')
  })
  
  try({
    ## selecting only the 10 genes which are DE in the most cells
    try(dev.off())
    pdf(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_percluster_DEgenes_subsetgenes_cluster_topgenes.pdf"),
        height = 4, width = 5)
    print(pheatmap::pheatmap(t(as((DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1][,order(apply(DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1], 2, function(i) sum(i > 0)), decreasing = T)[1:min(10, ncol(DE_Chrimsonposcells_lights_per_cluster_topgenesdcast)-1)]]), 'matrix'))))
    dev.off()
  })
  
  # print('here 5')
  DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes <- DE_Chrimsonposcells_lights_per_cluster_topgenes[DE_Chrimsonposcells_lights_per_cluster_topgenes$gene %in% names(which(table(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene) > 2)),]
  DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes$gene <- factor(DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes$gene, levels=names(which(table(DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes$gene) > 0)))
  table(DE_Chrimsonposcells_lights_per_cluster_topgenes$variable)
  
  cat('Plotting LFC of most important DEGs (1)...')
  # ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes %>% dplyr::filter(gene %in% c('sr', 'Hr38', 'CG14186', 'cbt', 'CG46385')), aes(y=gene, x=value))+geom_ridgeline()
  ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes, aes(y=gene, x=value))+geom_density_ridges()
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_percluster_DEgenes_subsetgenes_logfoldchange.pdf"), height = 3, width = 4)
  ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes,
         aes(y=gene, x=value))+geom_density_ridges()+
    geom_label(data = melt(table(DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes$gene)),
               aes(x=max(DE_Chrimsonposcells_lights_per_cluster_topgenes$value)-.5, y=Var1, label=paste0('n=',value)), hjust = 0,
               label.size = 0)
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_percluster_DEgenes_subsetgenes_logfoldchange_2.pdf"), height = 3, width = 4)
  
  cat('Plotting summary of most important DEGs (1)...')
  ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes,
         aes(y=gene, x=value))+geom_density_ridges()+
    geom_vline(xintercept = 0, lty='dashed', col='red')+
    geom_label(data = melt(table(DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes$gene)),
               aes(x=max(DE_Chrimsonposcells_lights_per_cluster_topgenes$value)-.5, y=Var1, label=paste0('n=',value)), hjust = 0,
               label.size = 0) #fill = NA
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_percluster_DEgenes_subsetgenes_logfoldchange_2_large.pdf"),
         height = min(0.7*length(levels(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene)), 6), width = 4)
  
  cat('Plotting summary of most important DEGs (2)...')
  ## summary of most important genes across clusters
  ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes[DE_Chrimsonposcells_lights_per_cluster_topgenes$gene %in% tail(levels(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene), n=10),], aes(y=gene, x=L1, fill=value))+
    # ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes, aes(y=gene, x=L1, fill=value))+
    scale_fill_gradientn(colours=c(bl,"white", re))+
    geom_tile()+ggtitle('Differentially expressed genes between lights ON/OFF per cell cluster\n(most shared genes)')+
    theme_bw()+labs(x='Cluster', y='Gene')+
    theme(axis.text.x = element_text(angle = 30, hjust=1))
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_percluster_DEgenes_subsetgenestop10.png"), height = 4, width = 8)
  
  # splitUMAPPlot(Chrimsonposcells, group.by='stim')
  # splitUMAPPlot(Chrimsonposcells, group.by='stim_light')
  
  
  cat('Computing number of DE genes...')
  ## num of DE genes per cluster
  dataset$num_light_DE_genes = table(DE_Chrimsonposcells_lights_per_cluster_topgenes$L1)[match(dataset$seurat_clusters, names(table(DE_Chrimsonposcells_lights_per_cluster_topgenes$L1)))]
  dataset$num_light_DE_genes_any_bool <- !is.na(dataset$num_light_DE_genes)
  dataset$num_light_DE_geneslogp1 <- log(dataset$num_light_DE_genes+1)
  table(dataset$num_light_DE_genes_any_bool)
  
  
  cat('Plotting features of interest...')
  ## UMAP of lights on/off coloured by feature
  splitDimRedPlot(dataset, group.by='stim_light', colour_by_feature='sr')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_sr.png"), height = 2.2, width = 4)
  splitDimRedPlot(dataset, group.by='stim_light', colour_by_feature='Hr38')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_Hr38.png"), height = 2.2, width = 4)
  splitDimRedPlot(dataset, group.by='stim_light', colour_by_feature='Hr38', dim_red = "tsne")
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_Hr38_TSNE.png"), height = 2.2, width = 4)
  splitDimRedPlot(dataset, group.by='stim_light', colour_by_feature='sr', dim_red = "tsne")
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_sr_TSNE.png"), height = 2.2, width = 4)
  
  cat('Plotting percentages of genes which are DE...')
  ## percentages of genes which are DE, in each cluster, including all cells, or separating by Chrimson cells
  dataset$num_light_DE_genes[is.na(dataset$num_light_DE_genes)] <- 0
  FeaturePlot(dataset, feature='num_light_DE_genes')+ggtitle('Number of DE genes')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_numDE_UMAP.png"), height = 2.8, width = 4)
  FeaturePlot(dataset, feature='num_light_DE_geneslogp1')+ggtitle('Number of DE genes')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_numDElogp1_UMAP.png"), height = 2.8, width = 4)
  FeaturePlot(dataset, feature='num_light_DE_genes_any_bool')+ggtitle('Presence of DE genes')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "lightONvsOFF_anyDE_UMAP.png"), height = 2.8, width = 4)
  
  # FeaturePlot(Chrimsonposcells, feature='num_light_DE_genes', reduction = 'pca')+ggtitle('Number of DE genes')
  # FeaturePlot(Chrimsonposcells, feature='num_light_DE_genes', reduction = 'tsne')+ggtitle('Number of DE genes')
  
  cat('Exploring power issues...')
  ## is there a power issue?
  num_DE_genes_cluster = cbind.data.frame(num_DE_genes=as.vector(table(factor(DE_Chrimsonposcells_lights_per_cluster_topgenes$L1, levels=sort(unique(dataset$seurat_clusters))))),
                                          cluster_size=as.vector(table(factor(dataset$seurat_clusters, levels=sort(unique(dataset$seurat_clusters))))),
                                          cluster=sort(unique(dataset$seurat_clusters)))

  ggplot(num_DE_genes_cluster, aes(x=cluster_size, y=num_DE_genes))+
    geom_point()+geom_smooth(method = "lm")+theme_bw()+labs(x='Size of cluster', y='Number of DE genes in cluster')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", "cor_numDEgenes_clustersize.png"), height = 3, width = 3)
  
  cat('Completed\n')
  # print('here 5')
  
}

python_like_select <- function(i, greplstr){ i[grepl(greplstr, i)] }

umap_plot_single_group <- function(seurat_obj, group.by, group_name, split_umap=F, na_as_other=F, ...){
  seurat_obj$cur <- seurat_obj@meta.data[,group.by] == group_name
  seurat_obj$cur[ seurat_obj$cur] <- group_name
  seurat_obj$cur[seurat_obj$cur == 'FALSE'] <- 'Other'
  if(na_as_other){
    seurat_obj$cur[is.na(seurat_obj$cur)] <- 'Other'
  }
  if(split_umap){
    splitDimRedPlot(seurat_obj, group.by='cur', ...)
  }else{
    UMAPPlot(seurat_obj, group.by='cur')
  }
}


transform_metadata_featureplot <- function(dataset, feature, trans='log2', ...){
  cat('Using umap by default')
  if(trans == 'log2')  dataset@meta.data[,'x'] <- log2(dataset@meta.data[,feature])
  else if(trans == 'log2p1')  dataset@meta.data[,'x'] <- log2(1+dataset@meta.data[,feature])
  FeaturePlot(dataset, features =  'x', ...)
}


transform_metadata_vlnplot <- function(dataset, feature, trans='log2', ...){
  if(trans == 'log2')  dataset@meta.data[,'x'] <- log2(dataset@meta.data[,feature])
  else if(trans == 'log2p1')  dataset@meta.data[,'x'] <- log2(1+dataset@meta.data[,feature])
  VlnPlot(dataset, features =  'x', ...)
}



Vlnplot_errorbar <- function(dataset, feature, normalised=F){
  if(normalised){
    x <- cbind.data.frame(counts=dataset@assays$RNA@data[feature,], ident=Idents(dataset))
  }else{
    x <- cbind.data.frame(counts=dataset@assays$RNA@counts[feature,], ident=Idents(dataset))
  }
  ggplot(x %>%
           group_by(ident) %>% summarise(median=median(counts), mean=mean(counts), sd=sd(counts, na.rm=T)) %>% arrange(median) %>% mutate(id = factor(ident, levels=rev(ident))),
         aes(x=id, y=median))+
    geom_point()+geom_line(aes(group=1))+
    geom_line(aes(y=mean, group=1), col='blue')+
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd,), col='blue')+
    theme(axis.text.x = element_text(angle = 30, hjust=1))+ggtitle(feature)+labs(x='Cluster', y='Median (black) and mean (blue) pm sd')
}

col_vector_2 <- c('#8ecae6','#219ebc','#023047','#ffb703','#fb8500', '#cdb4db','#ffc8dd','#ffafcc','#bde0fe','#a2d2ff',
                  '#75c08c', '#b6dc92','#d9f1d4', "#e63946","#f1faee","#a8dadc","#457b9d","#1d3557",
                  "#f4f1de","#e07a5f","#3d405b","#81b29a","#f2cc8f",
                  "#f72585","#7209b7","#3a0ca3","#4361ee","#4cc9f0",
                  "#ffadad","#ffd6a5","#fdffb6","#caffbf","#9bf6ff","#a0c4ff","#bdb2ff","#ffc6ff","#fffffc",
                  "#ffd6ff","#e7c6ff","#c8b6ff","#b8c0ff","#bbd0ff",
                  "#fbf8cc","#fde4cf","#ffcfd2","#f1c0e8","#cfbaf0","#a3c4f3","#90dbf4","#8eecf5","#98f5e1","#b9fbc0")

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

add_average_median_SD <- function(i, colname_mean_1='mean.Chrimson positive', colname_mean_2='mean.Chrimson negative'){
  a <- cbind(median=dcast(i, gene~L1, value.var = 'med'), mean=dcast(i, gene~L1, value.var = 'mean')[,-1],
             sd=dcast(i, gene~L1, value.var = 'sd')[,-1])
  a$col = as.factor(as.numeric(a[,colname_mean_1,] > a[,colname_mean_2]))
  a
}

plot_DE_top_cluster <- function(molten_df, top_n, ridge_plot=F){
  molten_df <- molten_df %>% group_by(L2, L1) %>% filter(row_number()%in%1:top_n)
  molten_df$L1L2 <- apply(molten_df[,c('L1', 'L2')],1, paste0, collapse='-')
  molten_df$signif <- molten_df$p_val_adj <= 0.05
  molten_df$L1 <- c('ChR-', 'ChR+')[match(molten_df$L1, c('chrimsonneg', 'chrimsonpos'))]
  
  if(ridge_plot){
    ggplot(molten_df,
           aes(y=factor(gene, levels=molten_df %>% 
                          group_by(gene) %>% 
                          dplyr::summarise(median(avg_log2FC)) %>% 
                          plyr::arrange(`median(avg_log2FC)`) %>% dplyr::select(gene) %>% unlist),
               x=avg_log2FC, fill=L1, alpha=0.2), lty=3)+geom_density_ridges()+theme_bw()+geom_vline(xintercept = 0)+
      scale_fill_jcolors(palette = "pal6")+labs(x='Average log-fold change across cluster', y='Gene', fill='')+
      theme(legend.position = "bottom")+guides(alpha="none")+
      facet_wrap(.~factor(signif, levels=c('TRUE', 'FALSE')), ncol=1)+ggtitle('LFC split by statistical differential abundance')
  }else{
    ggplot(molten_df, aes(y=factor(gene,levels=rev(unique(molten_df$gene))),
                          x=factor(L1L2, levels=rev(gtools::mixedsort(unique(L1L2)))), fill=avg_log2FC))+
      geom_tile()+
      geom_label(data = molten_df, aes(y=gene, x=L1L2, label=ifelse(avg_log2FC>0, '+', '-')), size=4, label.size = NA)+
      geom_point(data = molten_df %>% filter(signif), aes(y=gene, x=L1L2), col='black', size=5, shape=0)+
      theme(axis.text.x = element_text(angle = 30, hjust=1))+
      labs(x='Cluster', y='Gene')+
      facet_grid (.~ L1, scales = "free_x", space = "free_x")+
      scale_fill_gradient2( low = ("#c17889"), mid = "#f5f9a0", high = ("#9fdfc8"))
  }
}

predict_with_SingleR_consensus <- function(dataset, colname_new_metadata='singleRannotation', consensus_dataset){
  cat('Check that your datasets is of class: SingleCellExperiment\n')
  require(SingleR)
  if(consensus_dataset=='Aerts'){
    cat('Using one single reference cell of average expression per group\n')
    reference <- readRDS("~/projects/general/fly_midbrain_atlases/Atlas/AerstsSeuratConsensusCells.RDS")
  }else if(consensus_dataset=='Thirst_broad'){
    cat('Using one single reference cell of average expression per group\n')
    reference <- readRDS("~/projects/general/fly_midbrain_atlases/Atlas/thirstseuratConsensusCells_broad.RDS")
  }else if(consensus_dataset=='Thirst_specific'){
    cat('Using one single reference cell of average expression per group\n')
    reference <- readRDS("~/projects/general/fly_midbrain_atlases/Atlas/thirstseuratConsensusCells_specific.RDS")
  }
  annotation <- SingleR(test = dataset,
          ref = reference,
          labels = colnames(reference))
  dataset[[colname_new_metadata]] <- annotation
  return(dataset)
}
