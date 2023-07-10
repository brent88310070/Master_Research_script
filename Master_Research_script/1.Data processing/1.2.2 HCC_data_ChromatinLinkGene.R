library(Signac)
library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

#### Load data ####
HCC <- readRDS(file = './Human_Cerebral_Cortex/Cyc_N1_WNN_noLink_subset.rds')

#### Cut the coordinates ####
DimPlot(HCC, reduction = "wnn.umap", group.by = "celltype")
atacUMAP <- as.data.frame(HCC@reductions$wnn.umap@cell.embeddings)
atacUMAP <- dplyr::filter(atacUMAP, wnnUMAP_1 > -10 & wnnUMAP_1 < 3 & wnnUMAP_2 < 10 & wnnUMAP_2 > -10)

#remove cell
toKeep <- rownames(atacUMAP)
HCC <- HCC[,colnames(HCC) %in% toKeep]
DimPlot(HCC, reduction = "wnn.umap", group.by = "celltype")

#### separate nIPC and GluN1 ####
subset <- subset(x = HCC, subset = celltype == "nIPC/GluN1")
DimPlot(subset, reduction = "wnn.umap", group.by = "celltype")
DefaultAssay(subset) <- 'RNA'
FeaturePlot(subset, features = c("EOMES"), order = TRUE,reduction = "wnn.umap") + 
  ggtitle("EOMES") & scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "Reds")))

DefaultAssay(subset) <- 'ATAC' #ATAC for wnn
reduction <- "lsi" #wnn.umap
dims <- 1:10
subset <- FindNeighbors(object = subset, reduction = reduction, dims = dims)
subset <- FindClusters(object = subset, algorithm = 1, resolution = 1.5)
DimPlot(subset, reduction='wnn.umap', label = TRUE)
Meta_df <- subset@meta.data

library(car)
Meta_df$seurat_clusters <- recode(Meta_df$seurat_clusters ,"c('2', '6', '7', '11')='nIPC'")
Meta_df$seurat_clusters <- recode(Meta_df$seurat_clusters ,"c('0', '1', '3', '4', '5',
                                  '8', '9', '10', '12', '13', '14', '15')='GluN1'")
subset@meta.data$celltype <- Meta_df$seurat_clusters
DimPlot(subset, reduction='wnn.umap', group.by = "celltype", label = TRUE)
subset@meta.data$seurat_clusters <- NULL

reannotation_df <- data.frame(row.names = rownames(subset@meta.data), celltype=subset@meta.data$celltype)
annotation_df <- data.frame(row.names = rownames(HCC@meta.data), celltype=HCC@meta.data$celltype)
levels(annotation_df$celltype) <- c(levels(annotation_df$celltype), 'nIPC', 'GluN1')
library(tibble)
annotation_df <- annotation_df %>%
    rownames_to_column() %>%
    rows_update(reannotation_df %>% rownames_to_column(), by = "rowname") %>%
    column_to_rownames()

HCC@meta.data$celltype <- annotation_df$celltype
DimPlot(HCC, reduction='wnn.umap', group.by = "celltype", label = TRUE)

#### Combine N4+N5 ####
# celltype <- HCC@meta.data$celltype
# levels(celltype) <- c(levels(celltype), 'GluN4/5')
# celltype <- replace(celltype, celltype=="GluN4", "GluN4/5")
# celltype <- replace(celltype, celltype=="GluN5", "GluN4/5")
# HCC@meta.data$celltype <- celltype

#### Annotation ####
# extract gene annotations from EnsDb
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style
genome(annotation) <- "hg38"
annotation <- renameSeqlevels(annotation, mapSeqlevels(seqlevels(annotation), "UCSC"))
Annotation(HCC[["ATAC"]]) <- annotation

#### Subset cell ####
nIPC_GluN1 <- subset(x = HCC, subset = celltype == "nIPC" | celltype == "GluN1")
Cyc_nIPC <- subset(x = HCC, subset = celltype == "Cyc.Prog" | celltype == "nIPC")

HCC <- Cyc_nIPC
#### Correlation between ATAC and RNA ####
startTime <- Sys.time()
DefaultAssay(HCC) <- "ATAC"
HCC <- RegionStats(HCC, genome = BSgenome.Hsapiens.UCSC.hg38)
#Link peak to genes
HCC <- LinkPeaks(
  object = HCC,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  distance = 5e+05,
)

#get atac count
GetATAC_Count <- function(object, peaks){
  return(unname(rowSums(FetchData(object = HCC, vars = c(peaks)))))
}

#create atac peak filter df
df_gene <- unique(HCC@assays$ATAC@links$gene)
ATAC_df <- data.frame(matrix(ncol = length(colnames(HCC@assays$RNA@counts)), nrow = 0))
colnames(ATAC_df) <- colnames(HCC@assays$RNA@counts)
peaks_list <- HCC@assays$ATAC@links$peak[which(HCC@assays$ATAC@links$gene == df_gene[1])]
ATAC_df[1, ] <- GetATAC_Count(HCC, peaks_list)
for (i in 2:length(df_gene)) { 
  peaks_list <- HCC@assays$ATAC@links$peak[which(HCC@assays$ATAC@links$gene == df_gene[i])]
  ATAC_df[nrow(ATAC_df) + 1,] = GetATAC_Count(HCC, peaks_list)
} 
rownames(ATAC_df) <- df_gene
endTime <- Sys.time()
print(endTime - startTime)

#### Save data ####
HCC[['ATAC_filter']] <- CreateAssayObject(counts=ATAC_df)

# ----> save subset to RDS----
filename <- './Human_Cerebral_Cortex/DORC/Cyc_nIPC_filter_DORC.rds'
saveRDS(subset(x = HCC),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")
