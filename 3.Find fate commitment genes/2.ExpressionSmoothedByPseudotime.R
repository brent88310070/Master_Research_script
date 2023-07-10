library(Signac)
library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(SeuratWrappers)
library(monocle3)

#### Read data ####
CID <- readRDS(file = './filename.rds')

#### Remove Medulla part cells ####
Medulla_subset <- subset(x = CID, subset = celltype == "Medulla")
DimPlot(Medulla_subset, reduction = "wnn.umap", group.by = "celltype")
UMAP <- as.data.frame(Medulla_subset@reductions$wnn.umap@cell.embeddings)
UMAP <- dplyr::filter(UMAP, wnnUMAP_1 < 5 & wnnUMAP_1 > -10)
#Get remove cell
toRemove <- rownames(UMAP)
CID <- CID[,!colnames(CID) %in% toRemove]
DimPlot(CID, reduction = "wnn.umap", group.by = "celltype", label = TRUE) + NoLegend() 
rm(Medulla_subset)

#### Subset ####
Medulla <- subset(x = CID, subset = celltype == "Medulla" | celltype == "TAC" | celltype == "ORS")
DimPlot(Medulla, reduction = "wnn.umap", group.by = "celltype")
Cortex <- subset(x = CID, subset = celltype == "Hair Shaft-cuticle.cortex" | celltype == "TAC" | celltype == "ORS")
DimPlot(Cortex, reduction = "wnn.umap", group.by = "celltype")
IRS <- subset(x = CID, subset = celltype == "IRS" | celltype == "TAC" | celltype == "ORS")
DimPlot(Cortex, reduction = "wnn.umap", group.by = "celltype")

#### WNN pseudotime ####
PseudoT <- function(object){
  DefaultAssay(object) <- 'ATAC_filter'
  object@reductions$UMAP <- object@reductions$wnn.umap
  object.cds <- as.cell_data_set(object)
  object.cds <- cluster_cells(cds = object.cds, reduction_method = "UMAP")
  object.cds <- learn_graph(object.cds, use_partition = TRUE, close_loop = FALSE)
  object.cds <- order_cells(object.cds, reduction_method = "UMAP")
  WNNpseudoT <- object.cds@principal_graph_aux$UMAP$pseudotime
  return(WNNpseudoT)
}
ALL_PseudoT <- PseudoT(CID)
Medulla_PseudoT <- PseudoT(Medulla)
Cortex_PseudoT <- PseudoT(Cortex)
IRS_PseudoT <- PseudoT(IRS)

#Whole cell pseudotime plot
DefaultAssay(CID) <- 'ATAC_filter'
CID@reductions$UMAP <- CID@reductions$wnn.umap
CID.cds <- as.cell_data_set(CID)
CID.cds <- cluster_cells(cds = CID.cds, reduction_method = "UMAP")
CID.cds <- learn_graph(CID.cds, use_partition = TRUE, close_loop = FALSE)
CID.cds <- order_cells(CID.cds, reduction_method = "UMAP")
plot_cells(
  cds = CID.cds,
  color_cells_by = "pseudotime",
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE)

#### Particular genes expression smoothed by pseudotime (Loess) ####
minMax <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}
TimeSmooth <- function(object, PseudoT, gene, modal){
  DefaultAssay(object) <- modal
  df <- data.frame(PseudoT=PseudoT, expression=FetchData(object = object, vars = c(gene))[[gene]])
  smooth50 <- predict(loess(expression ~ PseudoT, data=df, span=.8))
  Nor_df <- data.frame(pseudoT=PseudoT, loess50=minMax(smooth50))
  return(Nor_df)
}

Smooth_plot <- function(object, PseudoT, TF, TG){
  TF_ATAC <- TimeSmooth(object, PseudoT, TF, "ATAC_filter")
  TF_RNA <- TimeSmooth(object, PseudoT, TF, "RNA")
  TG_ATAC <- TimeSmooth(object, PseudoT, TG, "ATAC_filter")
  TG_RNA <- TimeSmooth(object, PseudoT, TG, "RNA")
  PseudoT_df <- data.frame(pseudoT=PseudoT, celltype=object@meta.data$celltype)
  
  layout(matrix(c(1, 3, 
                  2, 4), nrow = 2, ncol = 2, byrow=TRUE),
         widths = c(6, 3), heights  = c(3, 1), respect = TRUE)
  
  par(mar = c(3, 3, 1.5, 0))
  plot(TF_ATAC$pseudoT, TF_ATAC$loess50, col = "red4",
       xlab = "", ylab = "", pch = 19, cex=.1) + 
      title(main=paste(TF, "→", TG), line=0.5)
  title(ylab="Normalized counts", line=2, cex.lab=1.2)
  title(xlab = "WNN Pseudotime", line=2, cex.lab=1.2)
  points(TF_RNA$pseudoT, TF_RNA$loess50, col='red', pch=19, cex=.1)
  points(TG_ATAC$pseudoT, TG_ATAC$loess50, col='blue4', pch=19, cex=.1)
  points(TG_RNA$pseudoT, TG_RNA$loess50, col='blue', pch=19, cex=.1)
  
  par(mar = c(1, 3, 0, 0))
  colors <- c("#FDAE61", "#D9EF8B", "#66BD63") 
  plot(PseudoT_df$pseudoT, c(rep(1, length(PseudoT_df$pseudoT))), main="",
       axes = FALSE, ylab = "", ylim = c(0.95, 1.05),
       col = colors[factor(PseudoT_df$celltype)], pch = "|", cex=.5)
  title(xlab = "Cell types", line=0, cex.lab=1.2)
  
  par(mar = c(0, 0, 0, 0))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend("left", legend=c(paste(TF, "ATAC"), paste(TF, "RNA"), 
         paste(TG, "ATAC"), paste(TG, "RNA")), bty = "n",
         cex=1, lwd=2, col=c('red4', 'red', 'blue4', 'blue'), x.intersp = 0.2, y.intersp = 0.8)

  par(mar = c(0, 0, 0, 0))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  cell_type <- levels(factor(PseudoT_df$celltype))
  legend("left", legend=c(levels(factor(PseudoT_df$celltype))), bty = "n",
         cex=1, lwd=2, col=colors, x.intersp = 0.2, y.intersp = 0.8,
         lty = c(NA, NA, NA), pch = c(15, 15, 15))
}

#Plot Normalization gene with ATAC & RNA
dev.off()
Smooth_plot(Cortex, Cortex_PseudoT, TF="Lef1", TG="Dsg4")


#### feature plot####
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
Gene_RNA_ATAC_plot('Gli2')


#### For Three plots ####
Smooth_plot_3_Genes <- function(object, PseudoT, g1, g2, g3){
  g1_ATAC <- TimeSmooth(object, PseudoT, g1, "ATAC_filter")
  g1_RNA <- TimeSmooth(object, PseudoT, g1, "RNA")
  g2_ATAC <- TimeSmooth(object, PseudoT, g2, "ATAC_filter")
  g2_RNA <- TimeSmooth(object, PseudoT, g2, "RNA")
  g3_ATAC <- TimeSmooth(object, PseudoT, g3, "ATAC_filter")
  g3_RNA <- TimeSmooth(object, PseudoT, g3, "RNA")
  PseudoT_df <- data.frame(pseudoT=PseudoT, celltype=object@meta.data$celltype)
  
  layout(matrix(c(1, 3, 
                  2, 4), nrow = 2, ncol = 2, byrow=TRUE),
         widths = c(6, 3), heights  = c(3, 1), respect = TRUE)
  
  par(mar = c(3, 3, 1.5, 0))
  plot(g1_ATAC$pseudoT, g1_ATAC$loess50, col = "red4",
       xlab = "", ylab = "", pch = 19, cex=.1) + 
    title(main=paste(g1, "→", g2, "→", g3), line=0.5)
  title(ylab="Normalized counts", line=2, cex.lab=1.2)
  title(xlab = "WNN Pseudotime", line=2, cex.lab=1.2)
  points(g1_RNA$pseudoT, g1_RNA$loess50, col='red', pch=19, cex=.1)
  points(g2_ATAC$pseudoT, g2_ATAC$loess50, col='blue4', pch=19, cex=.1)
  points(g2_RNA$pseudoT, g2_RNA$loess50, col='blue', pch=19, cex=.1)
  points(g3_ATAC$pseudoT, g3_ATAC$loess50, col='green4', pch=19, cex=.1)
  points(g3_RNA$pseudoT, g3_RNA$loess50, col='green', pch=19, cex=.1)
  
  par(mar = c(1, 3, 0, 0))
  colors <- c("#FDAE61", "#D9EF8B", "#66BD63") 
  plot(PseudoT_df$pseudoT, c(rep(1, length(PseudoT_df$pseudoT))), main="",
       axes = FALSE, ylab = "", ylim = c(0.95, 1.05),
       col = colors[factor(PseudoT_df$celltype)], pch = "|", cex=.5)
  title(xlab = "Cell types", line=0, cex.lab=1.2)
  
  par(mar = c(0, 0, 0, 0))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend("left", legend=c(paste(g1, "ATAC"), paste(g1, "RNA"), 
                          paste(g2, "ATAC"), paste(g2, "RNA"),
                          paste(g3, "ATAC"), paste(g3, "RNA")), bty = "n",
         cex=1, lwd=2, col=c('red4', 'red', 'blue4', 'blue', 'green4', 'green'), x.intersp = 0.2, y.intersp = 0.8)
  
  par(mar = c(0, 0, 0, 0))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  cell_type <- levels(factor(PseudoT_df$celltype))
  legend("left", legend=c(levels(factor(PseudoT_df$celltype))), bty = "n",
         cex=1, lwd=2, col=colors, x.intersp = 0.2, y.intersp = 0.8,
         lty = c(NA, NA, NA), pch = c(15, 15, 15))
}
Smooth_plot_3_Genes(Cortex, Cortex_PseudoT, g1="Lef1", g2="Hoxc13", g3="Wnt3")


