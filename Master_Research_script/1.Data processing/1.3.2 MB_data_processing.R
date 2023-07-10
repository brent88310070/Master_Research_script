library(Signac)
library(Seurat)
library(ggplot2)
library(magrittr)
library(data.table)
library(Matrix)
library(RColorBrewer)
library(patchwork)
library(EnsDb.Mmusculus.v79)

#### Read data ####
pathway <- './E18_Mouse_Brain/'
brain <- Read10X_h5(paste(pathway, "e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5", sep = ""), 
                  use.names = TRUE, unique.features = TRUE)

MB <- CreateSeuratObject(counts = brain$`Gene Expression`)
MB[['ATAC']] <- CreateChromatinAssay(counts = brain$Peaks, sep = c(":", "-"))
rm(brain)
gc()

#----> Annotation ----
# extract gene annotations from EnsDb
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style
genome(annotation) <- "mm10"
annotation <- renameSeqlevels(annotation, mapSeqlevels(seqlevels(annotation), "UCSC"))
Annotation(MB[["ATAC"]]) <- annotation

#### Read cell type file ####
celltype <- read.csv(file = paste(pathway, "celltype.csv", sep = ""))
celltype[celltype == "Astrocyte, Excitatory, Inhibitory"] <- "Progenitors"
celltype[celltype == "Excitatory, Inhibitory"] <- "Progenitors"

# Remove low quality cells
rm_barcodes <- setdiff(colnames(MB), celltype$X)
MB <- MB[,!colnames(MB) %in% rm_barcodes]

# Cell type annotation
MB@meta.data <- cbind(MB@meta.data, celltype$tree_states)
colnames(MB@meta.data)[which(names(MB@meta.data) == "celltype$tree_states")] <- "celltype"

# Cell type summary df
celltype.summary <- data.frame(summary(factor(MB@meta.data[[7]])))
colnames(celltype.summary)[1] <- "counts"
celltype.summary <- cbind(Celltype = rownames(celltype.summary), celltype.summary)
celltype.summary <- celltype.summary[order(-celltype.summary$counts),]
rownames(celltype.summary) <- 1:nrow(celltype.summary)

# Draw Violin Plot
VlnPlot(
  object = MB,
  features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC"),
  ncol = 4,
  pt.size = 0.01
)

#### Read UMAP file ####
umap_df <- read.csv(file = paste(pathway, "umap.csv", sep = ""))
umap_df$X <- NULL
rownames(umap_df) <- colnames(MB)
colnames(umap_df) <- c("umap_1", "umap_2")

umap_object <- CreateDimReducObject(embeddings = data.matrix(umap_df), key = "UMAP_", assay = "RNA", global = TRUE)
MB[['UMAP']] <- umap_object
DimPlot(MB, reduction = 'UMAP', label = TRUE, group.by = 'celltype') + NoLegend() + ggtitle("UMAP")

#### Save or Read file ####
#filename <- paste(pathway, "e18_mouse_brain_v2.rds", sep = "")
#saveRDS(MB ,file = filename)

MB <- readRDS(file = './E18_Mouse_Brain/e18_mouse_brain_v2.rds')

#### Normalization ####
#----> RNA Normalization ----
# Gene expression data processing
DefaultAssay(MB) <- "RNA"
MB <- FindVariableFeatures(MB, nfeatures = 3000)
MB <- NormalizeData(MB, normalization.method = "LogNormalize", scale.factor = 10000)
MB <- ScaleData(MB)
MB <- RunPCA(MB, npcs = 30)

#----> ATAC Normalization ----
DefaultAssay(MB) <- 'ATAC'
MB <- FindTopFeatures(MB, min.cutoff = 'q5')
MB <- RunTFIDF(MB)
MB <- RunSVD(MB)

#WNN UMAP gene (RNA-seq) 
Gene_RNA_ATAC_plot <- function(gene){
  DefaultAssay(MB) <- 'RNA'
  rna <- paste(gene, sep=" ")
  p1 <- FeaturePlot(MB, features = c(gene), order = TRUE, reduction = "UMAP") + ggtitle(rna) & scale_colour_gradientn(colours = (brewer.pal(n = 11, name = "Reds")))
  return(p1)
}
#Progenitors
Gene_RNA_ATAC_plot('Pax6')
#Astrocytes
Gene_RNA_ATAC_plot('Aqp4')
#Inhibitory neurons
Gene_RNA_ATAC_plot('Gad1')
#Excitatory neurons
Gene_RNA_ATAC_plot('Slc17a7')


#### Save Data ####
# ----> save subset to RDS----
# Progenitors and Excitatory
filename <- paste(pathway ,'Progenitors_Excitatory_noCorr_subset.rds', sep = "")
saveRDS(subset(x = MB, subset = celltype == "Progenitors" | celltype == "Excitatory" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# Progenitors and Inhibitory
filename <- paste(pathway ,'Progenitors_Inhibitory_noCorr_subset.rds', sep = "")
saveRDS(subset(x = MB, subset = celltype == "Progenitors" | celltype == "Inhibitory" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# Progenitors and Astrocyte
filename <- paste(pathway ,'Progenitors_Astrocyte_noCorr_subset.rds', sep = "")
saveRDS(subset(x = MB, subset = celltype == "Progenitors" | celltype == "Astrocyte" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# ALL
filename <- paste(pathway ,'ALL_noCorr_subset.rds', sep = "")
saveRDS(subset(x = MB, subset = celltype == "Progenitors" | celltype == "Astrocyte"| 
                 celltype == "Excitatory"| celltype == "Inhibitory" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")
