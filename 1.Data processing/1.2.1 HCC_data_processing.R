library(Signac)
library(Seurat)
library(ggplot2)
library(magrittr)
library(data.table)
library(Matrix)
library(EnsDb.Hsapiens.v86)
#set.seed(1234)

#### Read data ####
rna <- Read10X("./Human_Cerebral_Cortex/rna/", gene.column = 1)
cortex <- CreateSeuratObject(counts = rna)

atac <- Read10X("./Human_Cerebral_Cortex/atac/", gene.column = 1)
cortex[['ATAC']] <- CreateChromatinAssay(counts = atac, sep = c(":", "-"))
rm(rna)
rm(atac)
gc()


#### celltype annotation ####
cell_meta_data <- data.frame(fread(file="./Human_Cerebral_Cortex/GSE162170_multiome_cell_metadata.txt.gz", header=TRUE))
cortex@meta.data <- cbind(cortex@meta.data, cell_meta_data[, 66])
colnames(cortex@meta.data)[which(names(cortex@meta.data) == "cell_meta_data[, 66]")] <- "meta_cluster"

cortex <- SetIdent(cortex, value = cortex@meta.data$meta_cluster)
levels(cortex)

meta_cluster_celltype <- data.frame(fread(file="./Human_Cerebral_Cortex/GSE162170_multiome_cluster_names.txt.gz", header=TRUE))
cluster.celltype <- c("GluN5", "IN1", "nIPC/GluN1", "IN2", "SP", "GluN2", "IN3", 
                     "RG", "GluN4", "GluN3", "Cyc.Prog", "mGPC/OPC", "EC/Peric.", "mGPC/OPC")
names(cluster.celltype) <- levels(cortex)
cortex <- RenameIdents(cortex, cluster.celltype)

cortex@meta.data <- cbind(cortex@meta.data, celltype=Idents(cortex))

# Cell type summary df
celltype.summary <- data.frame(summary(factor(cortex@meta.data[[7]])))
colnames(celltype.summary)[1] <- "counts"
celltype.summary <- cbind(Celltype = rownames(celltype.summary), celltype.summary)
celltype.summary <- celltype.summary[order(-celltype.summary$counts),]
rownames(celltype.summary) <- 1:nrow(celltype.summary)

#----> Annotation ----
# extract gene annotations from EnsDb
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style
genome(annotation) <- "hg38"
annotation <- renameSeqlevels(annotation, mapSeqlevels(seqlevels(annotation), "UCSC"))
#seqlevelsStyle(annotation) <- "UCSC"
Annotation(cortex[["ATAC"]]) <- annotation

#### QC ####
Idents(cortex) <- 'orig.ident'
DefaultAssay(cortex) <- "RNA"
cortex[["percent.mt"]] <- PercentageFeatureSet(cortex, pattern = "^MT-")

DefaultAssay(cortex) <- "ATAC"
cortex$blacklist_fraction <- FractionCountsInRegion(
  object = cortex,
  assay = 'ATAC',
  regions = blacklist_hg38
)

# Draw Violin Plot
dev.off()
VlnPlot(
  object = cortex, features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC"),
  ncol = 4, pt.size = 0.01
)
VlnPlot(
  object = cortex, features = c("percent.mt", "blacklist_fraction"),
  ncol = 4, pt.size = 0.01
)

# Filter
#----> RNA ATAC QC ----
DefaultAssay(cortex) <- "RNA"

cortex <- subset(
  x = cortex,
  subset = blacklist_fraction < 0.03 &
    nCount_RNA > 500 &
    nCount_ATAC > 500
)

# Save or Read share-seq file
filename <- './Human_Cerebral_Cortex/Human_Cerebral_Cortex.rds'
saveRDS(cortex ,file = filename)

cortex <- readRDS(file = './Human_Cerebral_Cortex/Human_Cerebral_Cortex.rds')


#### Normalization ####
#----> RNA Normalization ----
# Gene expression data processing
DefaultAssay(cortex) <- "RNA"

cortex <- FindVariableFeatures(cortex, nfeatures = 2000)
cortex <- NormalizeData(cortex, normalization.method = "LogNormalize", scale.factor = 10000) #scale.factor = mean(share$nCount_RNA)

# RNA UMAP
cortex <- ScaleData(cortex)
cortex <- RunPCA(cortex, npcs = 50)

cortex <- RunUMAP(cortex, dims = 1:30, reduction.name = "umap.rna", min.dist = 0.8, n.neighbors = 50, metric = "cosine")
DimPlot(cortex, reduction = 'umap.rna', label = TRUE, group.by = 'celltype') + ggtitle("RNA UMAP") + NoLegend()

#----> ATAC Normalization ----
DefaultAssay(cortex) <- 'ATAC'
cortex <- FindTopFeatures(cortex, min.cutoff = 'q89.3')
cortex <- RunTFIDF(cortex)
cortex <- RunSVD(cortex)
DepthCor(cortex) + ggtitle("Human Cerebral Cortex")
cortex <- RunUMAP(cortex, metric = "cosine", dims = 2:10, reduction.name = 'umap.atac', min.dist = 0.6, n.neighbors = 50)
DimPlot(cortex, reduction = 'umap.atac', label = TRUE, group.by = 'celltype') + ggtitle("ATAC UMAP") + NoLegend()

# WNN
cortex <- FindMultiModalNeighbors(cortex, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 1:10))
cortex <- RunUMAP(cortex,  nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", 
                  n.neighbors = 50, metric = "cosine")
DimPlot(cortex, reduction = 'wnn.umap', label = TRUE, group.by = 'celltype') + ggtitle("WNN") + NoLegend()

p1 <- DimPlot(share, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size =4, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(share, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(share, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p1 &NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p2 &NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p3 &NoLegend() & theme(plot.title = element_text(hjust = 0.5))

# ----> save subset to RDS----
# nIPC/GluN1 and GluN2
filename <- './Human_Cerebral_Cortex/N1_N2_WNN_noLink_subset.rds'
saveRDS(subset(x = cortex, subset = celltype == "nIPC/GluN1" | celltype == "GluN2" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# GluN4 and GluN5
filename <- './Human_Cerebral_Cortex/N4_N5_WNN_noLink_subset.rds'
saveRDS(subset(x = cortex, subset = celltype == "GluN4" | celltype == "GluN5" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# GluN2, GluN4 and GluN5
filename <- './Human_Cerebral_Cortex/N2_N45_WNN_noLink_subset.rds'
saveRDS(subset(x = cortex, subset = celltype == "GluN2" | celltype == "GluN4" | celltype == "GluN5" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# GluN3, GluN4 and GluN5
filename <- './Human_Cerebral_Cortex/N3_N45_WNN_noLink_subset.rds'
saveRDS(subset(x = cortex, subset = celltype == "GluN3" | celltype == "GluN4" | celltype == "GluN5" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# cyc, Glu1
filename <- './Human_Cerebral_Cortex/Cyc_N1_WNN_noLink_subset.rds'
saveRDS(subset(x = cortex, subset = celltype == "Cyc.Prog" | celltype == "nIPC/GluN1" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# Glu2, Glu3
filename <- './Human_Cerebral_Cortex/N2_N3_WNN_noLink_subset.rds'
saveRDS(subset(x = cortex, subset = celltype == "GluN2" | celltype == "GluN3" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# Glu3, Glu4
filename <- './Human_Cerebral_Cortex/N3_N4_WNN_noLink_subset.rds'
saveRDS(subset(x = cortex, subset = celltype == "GluN3" | celltype == "GluN4" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")
