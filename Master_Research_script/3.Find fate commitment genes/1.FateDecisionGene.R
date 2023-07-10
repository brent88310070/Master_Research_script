library(ggvenn)
library(igraph)

#CuticleCortex Medulla IRS
CellType_On_1 <- readRDS("./share_seq/FateGenesNoTrans/CuticleCortex_On.Rda")
CellType_Off_1 <- readRDS("./share_seq/FateGenesNoTrans/CuticleCortex_Off.Rda")

CellType_On_2 <- readRDS("./share_seq/FateGenesNoTrans/Medulla_On.Rda")
CellType_Off_2 <- readRDS("./share_seq/FateGenesNoTrans/Medulla_Off.Rda")

threshold <- 0.05
CellType_On_1_filter_gene <- c(CellType_On_1[which(CellType_On_1$entropy_diff < -threshold),]$Gene)
CellType_Off_1_filter_gene <- c(CellType_Off_1[which(CellType_Off_1$entropy_diff > threshold),]$Gene)

CellType_On_2_filter_gene <- c(CellType_On_2[which(CellType_On_2$entropy_diff < -threshold),]$Gene)
CellType_Off_2_filter_gene <- c(CellType_Off_2[which(CellType_Off_2$entropy_diff > threshold),]$Gene)


# Mix (On + Off)
CellType_1_filter_gene <- c(CellType_On_1_filter_gene, CellType_Off_1_filter_gene)
CellType_2_filter_gene <- c(CellType_On_2_filter_gene, CellType_Off_2_filter_gene)

#### venn plot ####
l <- list("Cortex"=CellType_1_filter_gene, "Medulla"=CellType_2_filter_gene)
ggvenn(l)

CellType_1_diff <- setdiff(CellType_1_filter_gene, CellType_2_filter_gene)
CellType_2_diff <- setdiff(CellType_2_filter_gene, CellType_1_filter_gene)


#### Read TF TargetGene data ####
TFTG_df <- read.table(file = './mouse_tf_targetGene.tsv', sep = '\t', header = FALSE)
colnames(TFTG_df) <- c('TF', 'Target', 'Regulation', 'PMID')
TF <- unique(TFTG_df$TF)
TG <- unique(TFTG_df$Target)

#### Network of genes ####
Edge_list <- data.frame("TF"=TFTG_df$TF, "Target"=TFTG_df$Target)
G <- graph_from_data_frame(Edge_list, directed = TRUE)

Set <- CellType_1_diff
for(tf in intersect(Set, TF)){
  for(tg in intersect(Set, TG)){
    if(are.connected(G, tf, tg) == TRUE){
      print(paste(tf, "->", tg))
    }
  }
}

Set <- CellType_2_diff
for(tf in intersect(Set, TF)){
  for(tg in intersect(Set, TG)){
    if(are.connected(G, tf, tg) == TRUE){
      print(paste(tf, "->", tg))
    }
  }
}


#### GO enrichment ####
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
gene <- c()
GO_results <- enrichGO(gene = gene, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
GO_df <- as.data.frame(GO_results)
plot(barplot(GO_results, showCategory = 10))