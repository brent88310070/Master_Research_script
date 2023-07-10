library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(magrittr)
set.seed(1234)

#### Load Data####
#---> load processed data matrices for each assay----
rna <- Read10X("./SHARE/SHARE_RNA-seq", gene.column = 1)
atac <- Read10X("./SHARE_ATAC-seq", gene.column = 1)
fragments <- "./SHARE/GSM4156597_skin.late.anagen.atac.fragments.bed.gz"

# create a Seurat object and add the assays
share <- CreateSeuratObject(counts = rna)

share[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  fragments = fragments,
)

remove(rna, atac)
gc()

#----> Annotation ----
# extract gene annotations from EnsDb
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# change to UCSC style
genome(annotation) <- "mm10"
annotation <- renameSeqlevels(annotation, mapSeqlevels(seqlevels(annotation), "UCSC"))
#seqlevelsStyle(annotation) <- "UCSC"
Annotation(share[["ATAC"]]) <- annotation

#----> Celltype annotation ----
celltype_df <- read.table("./SHARE/GSM4156597_skin_celltype.txt", header = TRUE, sep = "\t")
colOrd <- rownames(FetchData(share,"ident"))
celltype_df <- celltype_df[match(noquote(colOrd), celltype_df$atac.bc),]

share@meta.data <- cbind(share@meta.data, celltype_df$celltype)
colnames(share@meta.data)[which(names(share@meta.data) == "celltype_df$celltype")] <- "celltype"

# Cell type summary df
celltype.summary <- data.frame(summary(factor(share@meta.data[[6]])))
colnames(celltype.summary)[1] <- "counts"
celltype.summary <- cbind(Celltype = rownames(celltype.summary), celltype.summary)
celltype.summary <- celltype.summary[order(-celltype.summary$counts),]
rownames(celltype.summary) <- 1:nrow(celltype.summary)

#### QC ####
DefaultAssay(share) <- "RNA"
share[["percent.mt"]] <- PercentageFeatureSet(share, pattern = "^mt-")

DefaultAssay(share) <- "ATAC"
share$blacklist_fraction <- FractionCountsInRegion(
  object = share,
  assay = 'ATAC',
  regions = blacklist_mm10
)

share <- TSSEnrichment(share)
share <- NucleosomeSignal(share)

# Filter
#----> RNA ATAC QC ----
DefaultAssay(share) <- "RNA"

share <- subset(
  x = share,
  subset = blacklist_fraction < 0.03 &
    nFeature_RNA > 100 &
    nFeature_RNA < 10000 &
    percent.mt < 2 &
    nCount_ATAC > 500 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2
)

# Draw Violin Plot
VlnPlot(
  object = share,
  features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "nFeature_ATAC"),
  ncol = 4,
  pt.size = 0.01
)

VlnPlot(
  object = share,
  features = c("percent.mt", "blacklist_fraction", 'nucleosome_signal', 'TSS.enrichment'),
  ncol = 4,
  pt.size = 0.01
)

# Remove mt genes
DefaultAssay(share) <- "RNA"
counts <- GetAssayData(share, assay = "RNA")
mt_gene <- rownames(counts[which(grepl("^mt-",rownames(counts))),])
counts <- counts[-(which(rownames(counts) %in% mt_gene)),]
dim(counts)
share@assays$RNA <- subset(x = share@assays$RNA, features = rownames(counts))

remove(counts)
gc()

# Save or Read share-seq file
#filename <- ''
#saveRDS(share ,file = filename)
share <- readRDS(file = './share_seq/mouse_skin.rds')

#### Normalization ####
#----> RNA Normalization ----
# Gene expression data processing
DefaultAssay(share) <- "RNA"

share <- FindVariableFeatures(share, nfeatures = 3000)
share <- NormalizeData(share, normalization.method = "LogNormalize", scale.factor = 10000) #scale.factor = mean(share$nCount_RNA)


# RNA UMAP
share <- ScaleData(share)
share <- RunPCA(share, npcs = 30)
share <- RunUMAP(share, dims = 1:30, reduction.name = "umap.rna")
DimPlot(share, reduction = 'umap.rna', label = TRUE, group.by = 'celltype') + NoLegend() + ggtitle("RNA UMAP")

#----> ATAC Normalization ----
DefaultAssay(share) <- 'ATAC'

share <- FindTopFeatures(share, min.cutoff = 'q5')
share <- RunTFIDF(share)
share <- RunSVD(share)
DepthCor(share) + ggtitle("SHARE-seq")
share <- RunUMAP(share, reduction = 'lsi', metric = "correlation", dims = 2:30, reduction.name = 'umap.atac')
DimPlot(share, reduction = 'umap.atac', label = TRUE, group.by = 'celltype') + NoLegend() + ggtitle("ATAC UMAP")

# WNN
share <- FindMultiModalNeighbors(share, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
share <- RunUMAP(share,  nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
DimPlot(share, reduction = 'wnn.umap', label = TRUE, group.by = 'celltype') + NoLegend() + ggtitle("WNN")

p1 <- DimPlot(share, reduction = "umap.rna", group.by = "celltype", label = TRUE, label.size =4, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(share, reduction = "umap.atac", group.by = "celltype", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(share, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 4, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p1 &NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p2 &NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p3 &NoLegend() & theme(plot.title = element_text(hjust = 0.5))

#### Peak gene link ####
library(BSgenome.Mmusculus.UCSC.mm10)

DefaultAssay(share) <- "ATAC"
share <- RegionStats(share, genome = BSgenome.Mmusculus.UCSC.mm10)
#Link peak to genes
share <- LinkPeaks(
  object = share,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  distance = 5e+05,
  n_sample = 100,
)
#### Gene activity score ####
DefaultAssay(share) <- 'ATAC'
gene.activities <- GeneActivity(share)

#----> To Binary data ----
gene.activities <- BinarizeCounts(gene.activities)
share[['ATAC_activity']] <- CreateAssayObject(counts = gene.activities)

# ----> save subset (with corr.) to RDS----
# Hair Follicles
filename <- './share_seq/Hair_Follicles_WNN_subset.rds'
saveRDS(subset(x = share, subset = celltype == "TAC" | celltype == "Hair Shaft-cuticle.cortex" | celltype == "Medulla"
               | celltype == "IRS" | celltype == "ORS"),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# TAC and Hair Shaft-cuticle.cortex
filename <- './share_seq/TAC_CuticleCortex_WNN_noCorr_subset.rds'
#SaveH5Seurat(subset(x = share, subset = celltype == "TAC" | celltype == "Hair Shaft-cuticle.cortex" ),filename = filename)
saveRDS(subset(x = share, subset = celltype == "TAC" | celltype == "Hair Shaft-cuticle.cortex" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# TAC and Medulla
filename <- './share_seq/TAC_Medulla_WNN_noCorr_subset.rds'
#SaveH5Seurat(subset(x = share, subset = celltype == "TAC" | celltype == "Medulla" ),filename = filename)
saveRDS(subset(x = share, subset = celltype == "TAC" | celltype == "Medulla" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# ----> save subset to RDS----
# TAC and Medulla
filename <- './share_seq/TAC_Medulla_WNN_subset.rds'
#SaveH5Seurat(subset(x = share, subset = celltype == "TAC" | celltype == "Medulla" ),filename = filename)
saveRDS(subset(x = share, subset = celltype == "TAC" | celltype == "Medulla" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# TAC and Cuticle&Cortex
filename <- './share_seq/TAC_CuticleCortex_WNN_subset.rds'
saveRDS(subset(x = share, subset = celltype == "TAC" | celltype == "Hair Shaft-cuticle.cortex" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# TAC and IRS
filename <- './share_seq/TAC_IRS_WNN_subset.rds'
saveRDS(subset(x = share, subset = celltype == "TAC" | celltype == "IRS" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")

# TAC and ORS
filename <- './share_seq/TAC_ORS_WNN_subset.rds'
saveRDS(subset(x = share, subset = celltype == "TAC" | celltype == "ORS" ),file = filename)
print(structure(file.size(filename), class = "object_size"), units = "Mb")
