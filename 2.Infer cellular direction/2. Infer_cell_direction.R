library(Signac)
library(Seurat)
library(ggplot2)
library(magrittr)
library(RColorBrewer)
library(patchwork)
library(dplyr)

#### File pathway ####
file_pathway <- "./share_seq/DORC/TAC_IRS_filter_WNN_DORC.rds"

##### Hyperparameter ####
Upstream <- "TAC"
Downstream <- "IRS"
percentile <- 0.25
sampling_percentage <- 0.6

#### Read data ####
ALL <- readRDS(file = file_pathway)

#### Filter by percentile ####
DefaultAssay(ALL) <- 'ATAC_filter'
remove_zero_quantile <- function(x){
  quantile(x[x!=0], c(percentile)) #percentile
}
binarize_by_threshold <- function(x){
  x > quantile_threshold
}
ATAC_count <- data.frame(ALL@assays$ATAC_filter@counts)
colnames(ATAC_count) <- colnames(ALL@assays$RNA)
quantile_threshold <- apply(ATAC_count, 1, remove_zero_quantile)
ATAC_count <- 1*apply(ATAC_count, 2, binarize_by_threshold)
ALL[['ATAC_Binarize_filter']] <- CreateAssayObject(counts = ATAC_count)

#colnames(ALL@meta.data)[6] <- "celltype"

#### Iteration START ####
START <- Sys.time()

times <- 100

#Hyperparameters
RNA_FC_value <- 1
RNA_p_value <- 0.05
OpenProbability <- 0.02
ATAC_FC_value <- 1.5

Result_df <- data.frame(matrix(nrow = times, ncol = 5))
for(TIME in 1:times){

  print(paste("##########", TIME, "##########"))
  #### Sampling Cells ####
  set.seed(TIME) #TIME
  toRemove <- sample(colnames(ALL), ceiling(dim(ALL)[2]*(1-sampling_percentage)))
  ALL_sampling <- ALL
  ALL_sampling <- ALL_sampling[, ! colnames(ALL_sampling) %in% toRemove]
  
  Upstream_subset <- subset(x = ALL_sampling, subset = celltype == Upstream)
  Downstream_subset <- subset(x = ALL_sampling, subset = celltype == Downstream)
  
  #### Get transition celltypes ####
  reduction <- "wnn.umap" #wnn.umap
  dims <- 1:2
  resolution <- 1.2

  #----> cluster larger celltype ----
  get_clusters_name <- function(large_object){
    larger_clusters_cellnames <- list()
    cluster_length <- length(levels(large_object@meta.data$seurat_clusters))
    for(i in 0:(cluster_length-1)){
      larger_clusters_cellnames[[i+1]] <- rownames(large_object@meta.data)[large_object@meta.data$seurat_clusters == i]
    }
    return(larger_clusters_cellnames)
  }
  if(dim(Upstream_subset)[2] >= dim(Downstream_subset)[2]){
    DefaultAssay(Upstream_subset) <- 'ATAC' #ATAC for wnn
    Upstream_subset <- FindNeighbors(object = Upstream_subset, reduction = reduction, dims = dims)
    Upstream_subset <- FindClusters(object = Upstream_subset, algorithm = 1, resolution = 1.2)
    larger_clusters_cellnames <- get_clusters_name(Upstream_subset)
  }else{
    DefaultAssay(Downstream_subset) <- 'ATAC' #ATAC for wnn
    Downstream_subset <- FindNeighbors(object = Downstream_subset, reduction = reduction, dims = dims)
    Downstream_subset <- FindClusters(object = Downstream_subset, algorithm = 1, resolution = 1.2)
    larger_clusters_cellnames <- get_clusters_name(Downstream_subset)
  }

  #----> Larger cell type: get each cluster center ----
  get_clusters_center <- function(large_object){
    cluster_length <- length(levels(large_object@meta.data$seurat_clusters))
    larger_clusters_center <- data.frame(matrix(nrow = cluster_length, ncol = length(dims)))
    cellsEmbeddings <- Embeddings(object = large_object[[reduction]])[, dims]
    for(i in 0:(cluster_length-1)){
      larger_clusters_center[i+1, ] <- colSums(cellsEmbeddings[larger_clusters_cellnames[[i+1]],]) / 
        length(larger_clusters_cellnames[[i+1]])
    }
    return(larger_clusters_center)
  }
  if(dim(Upstream_subset)[2] >= dim(Downstream_subset)[2]){
    larger_clusters_center <- get_clusters_center(Upstream_subset)
  }else{
    larger_clusters_center <- get_clusters_center(Downstream_subset)
  }

  #----> Get transition part center ----
  DefaultAssay(ALL_sampling) <- 'ATAC' #ATAC for wnn
  ALL_sampling <- FindNeighbors(object = ALL_sampling, reduction = reduction, dims = dims)
  ALL_sampling <- FindClusters(object = ALL_sampling, algorithm = 1, resolution = 1.2)
  
  cluster_celltype <- data.frame(Idents(ALL_sampling), ALL_sampling$celltype)
  cluster_df <- cluster_celltype[which(cluster_celltype$Idents.ALL_sampling.==0),]
  c1_num <- dim(cluster_df[which(cluster_df$ALL_sampling.celltype==Upstream),])[1]
  c2_num <- dim(cluster_df[which(cluster_df$ALL_sampling.celltype==Downstream),])[1]
  cluster_proportion_df <- data.frame(0, c1_num, c2_num, 
                                      round(c1_num/(c1_num+c2_num),2), round(c2_num/(c1_num+c2_num),2))
  colnames(cluster_proportion_df) <- c("Cluster", Upstream, Downstream, "Upstream_p", "Downstream_p")
  l <- length(unique(cluster_celltype$Idents.ALL_sampling.))-1
  for (i in 1:l){
    cluster_df <- cluster_celltype[which(cluster_celltype$Idents.ALL_sampling.==i),]
    c1_num <- dim(cluster_df[which(cluster_df$ALL_sampling.celltype==Upstream),])[1]
    c2_num <- dim(cluster_df[which(cluster_df$ALL_sampling.celltype==Downstream),])[1]
    cluster_proportion_df[nrow(cluster_proportion_df) + 1,] = c(i, c1_num, c2_num, round(c1_num/(c1_num+c2_num),2),
                                                                round(c2_num/(c1_num+c2_num),2))
  }
  transtion_cluster <- cluster_proportion_df[which.min(abs(0.5-cluster_proportion_df$Upstream_p)),]$Cluster
  clusters_cellnames <- rownames(ALL_sampling@meta.data)[ALL_sampling@meta.data$seurat_clusters == transtion_cluster]
  cellsEmbeddings <- Embeddings(object = ALL_sampling[[reduction]])[, dims]
  transition_cluster_center <- colSums(cellsEmbeddings[clusters_cellnames,]) / length(clusters_cellnames)
  
  
  #Find nearest cluster as transition
  calc_dis <- function(x){
    sqrt(sum((transition_cluster_center - x)^2))
  }
  dis <- apply(larger_clusters_center, 1, calc_dis)
  dis_df <- data.frame(Cluster=c(0:(length(dis)-1)), 
                       Cell_numbers=lengths(larger_clusters_cellnames), Distance=dis)
  dis_df <- dis_df[order(dis_df$Distance, decreasing = FALSE),]
  if(dim(Upstream_subset)[2] >= dim(Downstream_subset)[2]){
    clu <- sum(cumsum(dis_df$Cell_numbers) < dim(Downstream_subset)[2]) + 1
    keep_cluster <- paste(dis_df[1:clu, 1])
    Upstream_subset <- subset(x = Upstream_subset, idents = keep_cluster)
  }else{
    clu <- sum(cumsum(dis_df$Cell_numbers) < dim(Upstream_subset)[2]) + 1
    keep_cluster <- paste(dis_df[1:clu, 1])
    Downstream_subset <- subset(x = Downstream_subset, idents = keep_cluster)
  }
  # cell_types to infer direction (CID)
  KeepCell <- c(colnames(Upstream_subset), colnames(Downstream_subset))
  CID <- ALL_sampling[, colnames(ALL_sampling) %in% KeepCell]
  
  Upstream_subset <- subset(x = CID, subset = celltype == Upstream)
  Downstream_subset <- subset(x = CID, subset = celltype == Downstream)

  #Iterative nearest neighbor ( 10 clusters )
  Inn <- function(mat, clsize = 10, method=c('random','maxd', 'mind')){
    clsize.rle = rle(as.numeric(cut(1:nrow(mat), ceiling(nrow(mat)/clsize))))
    clsize = clsize.rle$lengths
    lab = rep(NA, nrow(mat))
    dmat = as.matrix(dist(mat))
    cpt = 1
    while(sum(is.na(lab)) > 0){
      lab.ii = which(is.na(lab))
      dmat.m = dmat[lab.ii,lab.ii]
      if(method[1]=='random'){
        ii = sample.int(nrow(dmat.m),1)
      } else if(method[1]=='maxd'){
        ii = which.max(rowSums(dmat.m))
      } else if(method[1]=='mind'){
        ii = which.min(rowSums(dmat.m))
      } else {
        stop('unknown method')
      }
      lab.m = rep(NA, length(lab.ii))
      lab.m[head(order(dmat.m[ii,]), clsize[cpt])] = cpt
      lab[lab.ii] = lab.m
      cpt = cpt + 1
    }
    if(any(is.na(lab))){
      lab[which(is.na(lab))] = cpt
    }
    return(lab)
  }
  
  #Upstream #wnn.umap umap.atac umap.rna
  DefaultAssay(Upstream_subset) <- 'ATAC' #ATAC for wnn
  if(dim(Upstream_subset@reductions$wnn.umap@cell.embeddings)[1] >= 100){
    Upstream_clsize <- ceiling(dim(Upstream_subset@reductions$wnn.umap@cell.embeddings)[1] / 10)
  }else{
    Upstream_clsize <- ceiling(dim(Upstream_subset@reductions$wnn.umap@cell.embeddings)[1] / 10) - 1
  }
  Upstream_Inn <- Inn(Upstream_subset@reductions$wnn.umap@cell.embeddings, clsize=Upstream_clsize, method='maxd')
  Upstream_subset@meta.data <- cbind(Upstream_subset@meta.data, Upstream_Inn)
  colnames(Upstream_subset@meta.data)[which(names(Upstream_subset@meta.data) == "Upstream_Inn")] <- "Inn10"
  #Downstream
  DefaultAssay(Downstream_subset) <- 'ATAC' #ATAC for wnn
  if(dim(Downstream_subset@reductions$wnn.umap@cell.embeddings)[1] >= 100){
    Downstream_clsize <- ceiling(dim(Downstream_subset@reductions$wnn.umap@cell.embeddings)[1] / 10)
  }else{
    Downstream_clsize <- ceiling(dim(Downstream_subset@reductions$wnn.umap@cell.embeddings)[1] / 10) - 1
  }
  Downstream_Inn <- Inn(Downstream_subset@reductions$wnn.umap@cell.embeddings, clsize=Downstream_clsize, method='maxd')
  Downstream_Inn <- Downstream_Inn + 10
  Downstream_subset@meta.data <- cbind(Downstream_subset@meta.data, Downstream_Inn)
  colnames(Downstream_subset@meta.data)[which(names(Downstream_subset@meta.data) == "Downstream_Inn")] <- "Inn10"
  
  Upstream_Inn_df <- data.frame(Inn10 = Upstream_subset@meta.data$Inn10, row.names = rownames(Upstream_subset@meta.data))
  Downstream_Inn_df <- data.frame(Inn10 = Downstream_subset@meta.data$Inn10, row.names = rownames(Downstream_subset@meta.data))
  Inn_df <- rbind(Upstream_Inn_df, Downstream_Inn_df)
  
  cellname <- rownames(FetchData(CID,"ident"))
  Inn_df <- Inn_df[match(noquote(cellname), rownames(Inn_df)),]
  CID@meta.data <- cbind(CID@meta.data, Inn_df)
  colnames(CID@meta.data)[which(names(CID@meta.data) == "Inn_df")] <- "Inn10"
  
  
  #Get sub cluster cell names
  cluster_cellnames <- list()
  for(i in 1:20){
    cluster_cellnames[[i]] <- rownames(CID@meta.data)[CID@meta.data$Inn10 == i]
  }
  subcluster_cell_number <- data.frame(cell_number = lengths(cluster_cellnames))
  
  #### Find DEGs ####
  DefaultAssay(CID) <- 'RNA'
  Idents(object = CID) <- CID@meta.data$celltype
  markers <- Seurat::FindAllMarkers(object = CID, test.use = "wilcox")
  markers <- markers[which(markers$avg_log2FC > log2(RNA_FC_value)),]
  markers <- markers[which(markers$p_val_adj < RNA_p_value),]
  
  #Remove some DEGs in df not in ATAC feature
  rm_DEGs <- setdiff(markers$gene, rownames(CID@assays$ATAC_Binarize_filter))
  markers <- markers[ ! markers$gene %in% rm_DEGs, ]
  
  #Create subset
  Turn_On_DEG <- markers[which(markers$cluster==Downstream),]
  Turn_Off_DEG <- markers[which(markers$cluster==Upstream),]
  length(Turn_On_DEG$gene)
  length(Turn_Off_DEG$gene)
  
  if (length(Turn_On_DEG$gene) == 0 | length(Turn_Off_DEG$gene) == 0) {
    next
  }
  
  ##### Filter gene with low atac-seq expression ####
  ATAC_filter_lowExpression <- function(Up_subset, Down_subset, genes){
    DefaultAssay(Up_subset) <- 'ATAC_Binarize_filter'
    Up_OpenCell <- colSums(FetchData(object = Up_subset, vars = c(genes)))
    DefaultAssay(Down_subset) <- 'ATAC_Binarize_filter'
    Down_OpenCell <- colSums(FetchData(object = Down_subset, vars = c(genes)))
    OpenProbability <- (Up_OpenCell + Down_OpenCell) / (dim(Up_subset)[2] + dim(Down_subset)[2])
    return(OpenProbability)
  }
  
  OpenProbability_vector <- ATAC_filter_lowExpression(Upstream_subset, Downstream_subset, markers$gene)
  DEG_ATAC_OpenProbability_df <- data.frame(Gene=markers$gene, OpenProbability=OpenProbability_vector)
  
  #### Filter gene with atac FC ####
  ATAC_FC <- function(situation, Up_subset, Down_subset, genes){
    DefaultAssay(Up_subset) <- 'ATAC_Binarize_filter'
    Up_percent <- colSums(FetchData(object = Up_subset, vars = c(genes))) / dim(Up_subset)[2]
    DefaultAssay(Down_subset) <- 'ATAC_Binarize_filter'
    Down_percent <- colSums(FetchData(object = Down_subset, vars = c(genes))) / dim(Down_subset)[2]
    if(situation == "On"){
      return(Down_percent/Up_percent)
    }
    if(situation == "Off"){
      return(Up_percent/Down_percent)
    }
  }
  #Turn on
  Turn_On_FC_vector <- ATAC_FC("On", Upstream_subset, Downstream_subset, Turn_On_DEG$gene)
  Turn_On_ATAC_FC_df <- data.frame(Gene=Turn_On_DEG$gene, FC=Turn_On_FC_vector)
  Turn_On_ATAC_FC_df <- Turn_On_ATAC_FC_df[is.finite(Turn_On_ATAC_FC_df$FC),]
  #Turn off
  Turn_Off_FC_vector <- ATAC_FC("Off", Upstream_subset, Downstream_subset, Turn_Off_DEG$gene)
  Turn_Off_ATAC_FC_df <- data.frame(Gene=Turn_Off_DEG$gene, FC=Turn_Off_FC_vector)
  Turn_Off_ATAC_FC_df <- Turn_Off_ATAC_FC_df[is.finite(Turn_Off_ATAC_FC_df$FC),]
  #ATAC DEG df
  DEG_ATAC_FC_df <- rbind(Turn_On_ATAC_FC_df, Turn_Off_ATAC_FC_df)
  DEG_ATAC_FC_df$FC[which(!is.finite(DEG_ATAC_FC_df$FC))] <- 0
  
  #### Filter low atac expression and low FC ####
  #Filter low open probability
  KeepGeneSet1 <- DEG_ATAC_OpenProbability_df[which(DEG_ATAC_OpenProbability_df$OpenProbability > OpenProbability),]$Gene
  #Filter low FC
  KeepGeneSet2 <- DEG_ATAC_FC_df[which(DEG_ATAC_FC_df$FC > ATAC_FC_value),]$Gene
  KeepGeneSet <- intersect(KeepGeneSet1, KeepGeneSet2)
  
  #Keep some DEGs in markers df
  RNA_markers <- markers
  markers <- markers[markers$gene %in% KeepGeneSet, ]
  length(markers$gene)
  
  #Create subset
  Turn_On_DEG <- markers[which(markers$cluster==Downstream),]
  Turn_Off_DEG <- markers[which(markers$cluster==Upstream),]
  
  #### Sampling Large subset ####
  cluster_KeepCell <- list()
  if (dim(Upstream_subset@assays$RNA@counts)[2] >= dim(Downstream_subset@assays$RNA@counts)[2]){
    sampling_size <- Downstream_clsize - 1
    for(i in 1:10){
      KeepCell_list <- sample(rownames(CID@meta.data[cluster_cellnames[[i]],]), sampling_size)
      cluster_KeepCell <- append(cluster_KeepCell, KeepCell_list)
    }
    Upstream_subset <- Upstream_subset[,colnames(Upstream_subset) %in% unlist(cluster_KeepCell)]
  }else{
    sampling_size <- Upstream_clsize - 1
    for(i in 11:20){
      KeepCell_list <- sample(rownames(CID@meta.data[cluster_cellnames[[i]],]), sampling_size)
      cluster_KeepCell <- append(cluster_KeepCell, KeepCell_list)
    }
    Downstream_subset <- Downstream_subset[,colnames(Downstream_subset) %in% unlist(cluster_KeepCell)]
  }
  
  #### Get sub cluster cell names ####
  #Upstream
  Upstream_cluster_cellnames <- list()
  for(i in 1:10){
    Upstream_cluster_cellnames[[i]] <- rownames(Upstream_subset@meta.data)[Upstream_subset@meta.data$Inn10 == i]
  }
  #Downstream
  Downstream_cluster_cellnames <- list()
  for(i in 11:20){
    Downstream_cluster_cellnames[[i]] <- rownames(Downstream_subset@meta.data)[Downstream_subset@meta.data$Inn10 == i]
  }
  
  #### Calculate each subcluster seq depth ####
  Calculate_subCluster_seqDepth <- function(upstream_object, downstream_object, whole_genes){
    #upstream
    DefaultAssay(upstream_object) <- 'ATAC_Binarize_filter'
    upstream_scale_factor <- c()
    for(i in 1:10){
      cluster_genes_df <- FetchData(upstream_object, vars = c(whole_genes), cells = Upstream_cluster_cellnames[[i]])
      cluster_genes_open_cell <- colSums(cluster_genes_df)
      cluster_mean_open_cell <- mean(cluster_genes_open_cell)
      upstream_scale_factor <- c(upstream_scale_factor, cluster_mean_open_cell)
    }
    #downstream
    DefaultAssay(downstream_object) <- 'ATAC_Binarize_filter'
    downstream_scale_factor <- c()
    for(i in 11:20){
      cluster_genes_df <- FetchData(downstream_object, vars = c(whole_genes), cells = Downstream_cluster_cellnames[[i]])
      cluster_genes_open_cell <- colSums(cluster_genes_df)
      cluster_mean_open_cell <- mean(cluster_genes_open_cell)
      downstream_scale_factor <- c(downstream_scale_factor, cluster_mean_open_cell)
    }
    return(round(c(upstream_scale_factor, downstream_scale_factor), digits = 3))
  }
  
  whole_genes <- intersect(rownames(CID@assays$RNA), rownames(CID@assays$ATAC_Binarize_filter))
  Subcluster_scaling_factor <- Calculate_subCluster_seqDepth(Upstream_subset, Downstream_subset, whole_genes)
  
  sd(Subcluster_scaling_factor[1:10]) #Upstream scaling factor SD
  sd(Subcluster_scaling_factor[11:20]) #Downstream scaling factor SD
  Scaling_factor_df <- data.frame(Scaling_factor=Subcluster_scaling_factor)
  
  #### Calculate Entropy ####
  #Entropy function
  shannon_entropy <- function(p){
    round(-sum(na.omit(p * log2(p))),3)
  }
  Calculate_entropy <- function(upstream_object, downstream_object, genes){
    #upstream
    DefaultAssay(upstream_object) <- 'ATAC_Binarize_filter'
    upstream_p_df <- data.frame(matrix(nrow = 10, ncol = length(genes)))
    for(i in 1:10){
      open_cell <- colSums(FetchData(upstream_object, vars = c(genes), 
                                     cells = Upstream_cluster_cellnames[[i]]))
      p <- open_cell / length(Upstream_cluster_cellnames[[i]])
      upstream_p_df[i, ] <- p
    }
    scale_upstream_p_df <- upstream_p_df / Subcluster_scaling_factor[1:10]
    Nor_upstream_p_df <- sweep(scale_upstream_p_df, 2, colSums(scale_upstream_p_df), "/")
    upstream_entropy <- apply(Nor_upstream_p_df, 2, shannon_entropy)
    
    #downstream
    DefaultAssay(downstream_object) <- 'ATAC_Binarize_filter'
    downstream_p_df <- data.frame(matrix(nrow = 10, ncol = length(genes)))
    for(i in 1:10){
      open_cell <- colSums(FetchData(downstream_object, vars = c(genes), 
                                     cells = Downstream_cluster_cellnames[[i+10]]))
      p <- open_cell / length(Downstream_cluster_cellnames[[i+10]])
      downstream_p_df[i, ] <- p
    }
    scale_downstream_p_df <- downstream_p_df / Subcluster_scaling_factor[11:20]
    Nor_downstream_p_df <- sweep(scale_downstream_p_df, 2, colSums(scale_downstream_p_df), "/")
    downstream_entropy <- apply(Nor_downstream_p_df, 2, shannon_entropy)
    
    return(data.frame(Gene=genes, Upstream_H=unname(upstream_entropy), Downstream_H=unname(downstream_entropy)))
  }
  
  #Calculate DEGs entropy and entropy difference 
  entropy_df <- Calculate_entropy(Upstream_subset, Downstream_subset, markers$gene)
  entropy_df['entropy_diff'] <- entropy_df$Upstream_H-entropy_df$Downstream_H
  
  #Turn on DEGs entropy
  Turn_on_entropy_df <- filter(entropy_df, Gene %in% Turn_On_DEG$gene)
  Turn_on_entropy_df <- Turn_on_entropy_df[rowSums(is.na(Turn_on_entropy_df)) == 0, ] 
  
  #Turn off DEGs entropy
  Turn_off_entropy_df <- filter(entropy_df, Gene %in% Turn_Off_DEG$gene)
  Turn_off_entropy_df <- Turn_off_entropy_df[rowSums(is.na(Turn_off_entropy_df)) == 0, ] 
  
  On_mean <- abs(mean(Turn_on_entropy_df$entropy_diff))
  Off_mean <- abs(mean(Turn_off_entropy_df$entropy_diff))
  On_num <- length(Turn_on_entropy_df$Gene)
  Off_num <- length(Turn_off_entropy_df$Gene)
  
  Result_df[TIME,] <- c(On_mean, Off_mean, On_num, Off_num, (On_mean > Off_mean))
}
colnames(Result_df) <- c("On_H", "Off_H", "On_num", "Off_num", "Ans.")

END <- Sys.time()
print(END - START)

# Print accuracy
print(paste("Accuracy: ", sum(Result_df$Ans.), "%"))

mean(Result_df$On_num)
mean(Result_df$Off_num)

sum(na.omit(Result_df$Ans.)) / length(na.omit(Result_df$Ans.))
mean(na.omit(Result_df$On_num))
mean(na.omit(Result_df$Off_num))

#WNN UMAP gene (RNA + ATAC-seq) 
Gene_RNA_ATAC_plot <- function(gene){
  DefaultAssay(CID) <- 'RNA'
  rna <- paste(gene, "(RNA)", sep=" ")
  p1 <- FeaturePlot(CID, features = c(gene), order = TRUE,reduction = "wnn.umap") + ggtitle(rna) & scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "Reds")))
  DefaultAssay(CID) <- 'ATAC_Binarize_filter'
  atac <- paste(gene, "(ATAC)", sep=" ")
  p2 <- FeaturePlot(CID, features = c(gene), order = TRUE,reduction = "wnn.umap") + ggtitle(atac) & scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "Reds"))) 
  p2 <- p2 + theme(legend.position="none")
  p3 <- p1 + p2 + plot_layout(ncol=2)
  return(p3)
}
Gene_RNA_ATAC_plot('WDR43')
