library(Signac)
library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(BSgenome.Mmusculus.UCSC.mm10)

#### Load data ####
pathway <- '/media/data/single_cell/brent10070/E18_Mouse_Brain/'
MB <- readRDS(file = paste(pathway ,'ALL_noCorr_subset.rds', sep = ''))

#### Cut the coordinates ####
DimPlot(MB, reduction = "UMAP", group.by = "celltype")
UMAP <- as.data.frame(MB@reductions$UMAP@cell.embeddings)
UMAP <- dplyr::filter(UMAP, UMAP_1 > 1 & UMAP_1 < 15 & UMAP_2 < 15 & UMAP_2 > 0)

#remove cell
toKeep <- rownames(UMAP)
MB <- MB[,colnames(MB) %in% toKeep]
DimPlot(MB, reduction = "UMAP", group.by = "celltype")

# rm_barcodes <- c('ACAGCCGGTACGGTAC-1', 'GGCTCACAGGAACCAA-1')
# MB <- MB[,!colnames(MB) %in% rm_barcodes]

#### Correlation between ATAC and RNA ####
startTime <- Sys.time()
DefaultAssay(MB) <- "ATAC"
MB <- RegionStats(MB, genome = BSgenome.Mmusculus.UCSC.mm10)
#Link peak to genes
MB <- LinkPeaks(
  object = MB,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  distance = 5e+05,
)

#get atac count
GetATAC_Count <- function(object, peaks){
  return(unname(rowSums(FetchData(object = MB, vars = c(peaks)))))
}

#create atac peak filter df
df_gene <- unique(MB@assays$ATAC@links$gene)
ATAC_df <- data.frame(matrix(ncol = length(colnames(MB@assays$RNA@counts)), nrow = 0))
colnames(ATAC_df) <- colnames(MB@assays$RNA@counts)
peaks_list <- MB@assays$ATAC@links$peak[which(MB@assays$ATAC@links$gene == df_gene[1])]
ATAC_df[1, ] <- GetATAC_Count(MB, peaks_list)
for (i in 2:length(df_gene)) { 
  peaks_list <- MB@assays$ATAC@links$peak[which(MB@assays$ATAC@links$gene == df_gene[i])]
  ATAC_df[nrow(ATAC_df) + 1,] = GetATAC_Count(MB, peaks_list)
} 
rownames(ATAC_df) <- df_gene
endTime <- Sys.time()
print(endTime - startTime)

#### Save data ####
MB[['ATAC_filter']] <- CreateAssayObject(counts=ATAC_df)

# ----> save subset to RDS----
filename <- paste(pathway ,'DORC/ALL_filter_DORC.rds', sep = "")
saveRDS(subset(x = MB),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")
