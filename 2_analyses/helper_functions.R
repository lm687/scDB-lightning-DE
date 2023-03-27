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

add_dimred <- function(seurat_obj, nPCs=60, resolution_clusters=4, compute_neighbours=T){
  cat('Scaling...\n')
  seurat_obj <- ScaleData(seurat_obj)
  cat('Running PCA...\n')
  seurat_obj <- RunPCA(seurat_obj, npcs = 60)
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
  cat('Here\n')
  res <- lapply(clusters, function(cluster_it){
    idx_select <- which(seurat_obj@meta.data[,cluster_name] == cluster_it)
    print(idx_select)
    subset_cells <- subset(seurat_obj, cells = idx_select)
    cat('Number of cells in each DE group: ', length(unique(subset_cells@meta.data[,group.by_arg])), '\n')
    if( length(unique(subset_cells@meta.data[,group.by_arg])) == 1 ){
      cat('There is perfect separation between DE groups in cluster', cluster_it, '\n')
      warning('There is perfect separation between DE groups in cluster ', cluster_it)
      ## there is perfect separation
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


give_DE_analysis_lightON_lightOFF <- function(dataset, dataset_name, name_clusters=NULL, add_name=''){
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
  DE_Chrimsonposcells_lights <- Seurat::FoldChange(dataset, ident.1='Lights on', ident.2='Lights off', group.by='stim_light')
  DE_Chrimsonposcells_lights
  DE_Chrimsonposcells_lights_topgenes <- give_top_logf_genes(DE_Chrimsonposcells_lights)
  xtable::xtable(DE_Chrimsonposcells_lights_topgenes)
  
  system(paste0("mkdir -p ", folder_results, input_objs, "/", dataset_name, add_name))
  
  FeaturePlot(dataset, features = rownames(DE_Chrimsonposcells_lights_topgenes))
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_featureplot.png"), height = 6, width = 8)
  
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
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes.png"),
         height = min(length(levels(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene)), 15), width = 8)
  
  
  ## to find clustering
  select_two_or_more_active <- function(i){
    i[,colSums(i>0) > 2]
  }
  
  try({
    try(dev.off())
    pdf(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes_subsetgenes_cluster.pdf"), height = 3, width = 4)
    # print(pheatmap::pheatmap(t(as(select_two_or_more_active(DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1]), 'matrix'))))
    print(pheatmap::pheatmap(t(as((DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1]), 'matrix'))))
    dev.off()
    # print('here 4')
  })
  
  try({
    try(dev.off())
    pdf(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes_subsetgenes_cluster_large.pdf"),
        height = min(length(levels(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene)), 15), width = 8)
    # print(pheatmap::pheatmap(t(as(select_two_or_more_active(DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1]), 'matrix'))))
    print(pheatmap::pheatmap(t(as((DE_Chrimsonposcells_lights_per_cluster_topgenesdcast[,-1]), 'matrix'))))
    dev.off()
    # print('here 4')
  })
  
  try({
    ## selecting only the 10 genes which are DE in the most cells
    try(dev.off())
    pdf(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes_subsetgenes_cluster_topgenes.pdf"),
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
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes_subsetgenes_logfoldchange.pdf"), height = 3, width = 4)
  ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes,
         aes(y=gene, x=value))+geom_density_ridges()+
    geom_label(data = melt(table(DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes$gene)),
               aes(x=max(DE_Chrimsonposcells_lights_per_cluster_topgenes$value)-.5, y=Var1, label=paste0('n=',value)), hjust = 0,
               label.size = 0)
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes_subsetgenes_logfoldchange_2.pdf"), height = 3, width = 4)
  
  cat('Plotting summary of most important DEGs (1)...')
  ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes,
         aes(y=gene, x=value))+geom_density_ridges()+
    geom_vline(xintercept = 0, lty='dashed', col='red')+
    geom_label(data = melt(table(DE_Chrimsonposcells_lights_per_cluster_topgenes_subset_genes$gene)),
               aes(x=max(DE_Chrimsonposcells_lights_per_cluster_topgenes$value)-.5, y=Var1, label=paste0('n=',value)), hjust = 0,
               label.size = 0) #fill = NA
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes_subsetgenes_logfoldchange_2_large.pdf"),
         height = min(0.7*length(levels(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene)), 6), width = 4)
  
  cat('Plotting summary of most important DEGs (2)...')
  ## summary of most important genes across clusters
  ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes[DE_Chrimsonposcells_lights_per_cluster_topgenes$gene %in% tail(levels(DE_Chrimsonposcells_lights_per_cluster_topgenes$gene), n=10),], aes(y=gene, x=L1, fill=value))+
    # ggplot(DE_Chrimsonposcells_lights_per_cluster_topgenes, aes(y=gene, x=L1, fill=value))+
    scale_fill_gradientn(colours=c(bl,"white", re))+
    geom_tile()+ggtitle('Differentially expressed genes between lights ON/OFF per cell cluster\n(most shared genes)')+
    theme_bw()+labs(x='Cluster', y='Gene')+
    theme(axis.text.x = element_text(angle = 30, hjust=1))
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_percluster_DEgenes_subsetgenestop10.png"), height = 4, width = 8)
  
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
  splitUMAPPlot(dataset, group.by='stim_light', colour_by_feature='sr')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_sr.png"), height = 2.2, width = 4)
  splitUMAPPlot(dataset, group.by='stim_light', colour_by_feature='Hr38')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_Hr38.png"), height = 2.2, width = 4)
  splitUMAPPlot(dataset, group.by='stim_light', colour_by_feature='Hr38', dim_red = "tsne")
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_Hr38_TSNE.png"), height = 2.2, width = 4)
  splitUMAPPlot(dataset, group.by='stim_light', colour_by_feature='sr', dim_red = "tsne")
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_sr_TSNE.png"), height = 2.2, width = 4)
  
  cat('Plotting percentages of genes which are DE...')
  ## percentages of genes which are DE, in each cluster, including all cells, or separating by Chrimson cells
  dataset$num_light_DE_genes[is.na(dataset$num_light_DE_genes)] <- 0
  FeaturePlot(dataset, feature='num_light_DE_genes')+ggtitle('Number of DE genes')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_numDE_UMAP.png"), height = 2.8, width = 4)
  FeaturePlot(dataset, feature='num_light_DE_geneslogp1')+ggtitle('Number of DE genes')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_numDElogp1_UMAP.png"), height = 2.8, width = 4)
  FeaturePlot(dataset, feature='num_light_DE_genes_any_bool')+ggtitle('Presence of DE genes')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_lightONvsOFF_anyDE_UMAP.png"), height = 2.8, width = 4)
  
  # FeaturePlot(Chrimsonposcells, feature='num_light_DE_genes', reduction = 'pca')+ggtitle('Number of DE genes')
  # FeaturePlot(Chrimsonposcells, feature='num_light_DE_genes', reduction = 'tsne')+ggtitle('Number of DE genes')
  
  cat('Exploring power issues...')
  ## is there a power issue?
  num_DE_genes_cluster = cbind.data.frame(num_DE_genes=as.vector(table(factor(DE_Chrimsonposcells_lights_per_cluster_topgenes$L1, levels=sort(unique(dataset$seurat_clusters))))),
                                          cluster_size=as.vector(table(factor(dataset$seurat_clusters, levels=sort(unique(dataset$seurat_clusters))))),
                                          cluster=sort(unique(dataset$seurat_clusters)))
  num_DE_genes_cluster
  
  ggplot(num_DE_genes_cluster, aes(x=cluster_size, y=num_DE_genes))+
    geom_point()+geom_smooth(method = "lm")+theme_bw()+labs(x='Size of cluster', y='Number of DE genes in cluster')
  ggsave(paste0(folder_results, input_objs, "/", dataset_name, add_name, "/", dataset_name, "_cor_numDEgenes_clustersize.png"), height = 3, width = 3)
  
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
