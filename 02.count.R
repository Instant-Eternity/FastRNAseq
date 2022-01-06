##########################################################################
# File Name: 02.count.R
# Author: Instant-Eternity
# mail: hunterfirstone@i.smu.edu.cn
# Created Time: Thu 03 Jun 2021 11:34:14 AM CST
#########################################################################
#!/usr/bin/env Rscript
# Clean work space
rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

getwd()
setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/19_LJX_RNAseq/08.fpkm")

Raw_count <- read.table("/bigdata/wangzhang_guest/chenpeng_project/06_result/19_LJX_RNAseq/all.count",header = T)

Raw_count_filt <- Raw_count

ENSEMBL <- gsub("\\.\\d*", "", Raw_count_filt$gene_id)
row.names(Raw_count_filt) <- ENSEMBL
head(Raw_count_filt)

dir.create("/bigdata/wangzhang_guest/chenpeng_project/06_result/19_LJX_RNAseq/09.MergeCount")
setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/19_LJX_RNAseq/09.MergeCount")

saveRDS(Raw_count_filt,"Raw_count_filt")
Raw_count_filt <- readRDS("Raw_count_filt")
write.csv(Raw_count_filt, file="Readcount_fpkm.csv")

library("org.Mm.eg.db")

g2s <- unique(toTable(org.Mm.egSYMBOL))
g2e <- unique(toTable(org.Mm.egENSEMBL))
s2e <- merge(g2e, g2s, by = 'gene_id')

table(row.names(Raw_count_filt) %in% s2e$ensembl)
Raw_count_filt <- Raw_count_filt[row.names(Raw_count_filt) %in% s2e$ensembl,]
Raw_count_filt$gene_symbol <- s2e[match(row.names(Raw_count_filt), s2e$ensembl),3]
length(unique(Raw_count_filt$gene_id))

saveRDS(Raw_count_filt,"Raw_count_with_symbol")
Raw_count_filt <- readRDS("Raw_count_with_symbol")
write.csv(Raw_count_filt, file="Readcount_fpkm_with_symbol.csv")

library(DESeq2)
Raw_count_filt <- read.csv("Readcount_fpkm_with_symbol.csv",row.names=1,header=T)

#database_count <- as.matrix(Raw_count_filt[,c(2,5,8,11,14,17,20,23,26,29,32,35,38,41,44)])
#database_fpkm <- as.matrix(Raw_count_filt[,c(3,6,9,12,15,18,21,24,27,30,33,36,39,42,45)])
#database_tpm <- as.matrix(Raw_count_filt[,c(4,7,10,13,16,19,22,25,28,31,34,37,40,43,46)])

Count_Deg <- function(Count,colData,condition,ref_group){
    dds <- DESeqDataSetFromMatrix(Count, colData, design = ~condition)
    dds <- dds[rowSums(counts(dds)) > 1,]
    dds$condition <- relevel(dds$condition, ref = ref_group)
    dds <- estimateSizeFactors(dds)
    dds <- DESeq(dds)
    res <- results(dds)
    resordered <- res[order(res$padj),]
    summary(res)
    return (res)
}

Count_shame_clp <- as.matrix(Raw_count_filt[,c(2,5,8,11,14,32,35,38,41,44)])
condition <- factor(c(rep('Shame',each=5),rep('Clp',each=5)))
colData <- data.frame(row.names = colnames(Count_shame_clp[1:10]),condition)
result_shame_clp <- Count_Deg(Count_shame_clp,colData,condition,'Shame')

Count_clp_HDCA <- as.matrix(Raw_count_filt[,c(17,20,23,26,29,32,35,38,41,44)])
condition <- factor(c(rep('HDCA',each=5),rep('Clp',each=5)))
colData <- data.frame(row.names = colnames(Count_clp_HDCA[1:10]),condition)
result_clp_HDCA <- Count_Deg(Count_clp_HDCA,colData,condition,'Clp')

dir.create("/bigdata/wangzhang_guest/chenpeng_project/06_result/19_LJX_RNAseq/10.DEG")
setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/19_LJX_RNAseq/10.DEG")

write.csv(result_shame_clp,file = "ALL_shame_vs_clp.csv")
table(result_shame_clp$padj<0.05)

write.csv(result_clp_HDCA,file = "ALL_clp_vs_HDCA.csv")
table(result_clp_HDCA$padj<0.05)

diff_gene_deseq2 <-subset(result_shame_clp,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
dim(diff_gene_deseq2)
head(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "DEG_shame_vs_clp.csv")

diff_gene_deseq2 <-subset(result_clp_HDCA,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
dim(diff_gene_deseq2)
head(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "DEG_clp_HDCA.csv")

library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)

setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/19_LJX_RNAseq/10.DEG")
#sig.gene <- read.csv("DEG_shame_vs_clp.csv")
sig.gene <- read.csv("DEG_clp_HDCA.csv")

dir.create("/bigdata/wangzhang_guest/chenpeng_project/06_result/19_LJX_RNAseq/11.GO")
setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/19_LJX_RNAseq/11.GO")
head(sig.gene)
#gene <- row.names(sig.gene)
gene <- sig.gene[,1]
gene.df<-bitr(gene, fromType = "ENSEMBL",
              toType = c("SYMBOL","ENTREZID"),
              OrgDb = org.Mm.eg.db)
display_number = c(15, 15, 15)
ego_MF <- enrichGO(OrgDb="org.Mm.eg.db",
                   gene = gene.df$ENSEMBL,
                   keyType = 'ENSEMBL',
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   ont = "MF",
                   pAdjustMethod = "BH",
                   readable=TRUE)
ego_result_MF <- as.data.frame(ego_MF)[1:display_number[1], ]
ego_CC <- enrichGO(OrgDb="org.Mm.eg.db",
                   gene = gene.df$ENSEMBL,
                   keyType = 'ENSEMBL',
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   ont = "CC",
                   pAdjustMethod = "BH",
                   readable=TRUE)
ego_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_BP <- enrichGO(OrgDb="org.Mm.eg.db",
                   gene = gene.df$ENSEMBL,
                   keyType = 'ENSEMBL',
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   readable=TRUE)
ego_result_BP <- na.omit(as.data.frame(ego_BP))[1:display_number[3], ]
go_enrich_df <- data.frame(ID=c(ego_result_BP$ID, ego_result_CC$ID, ego_result_MF$ID),
                           Description=c(ego_result_BP$Description, ego_result_CC$Description, ego_result_MF$Description),
                           GeneNumber=c(ego_result_BP$Count, ego_result_CC$Count, ego_result_MF$Count),
                           type=factor(c(rep("biological process", display_number[1]), rep("cellular component", display_number[2]),rep("molecular function", display_number[3])), levels=c("molecular function", "cellular component", "biological process")))

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))

## numbers as data on x axis
go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))
## shorten the names of GO terms
shorten_names <- function(x, n_word=4, n_char=45){
    if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 45))
        {
            if (nchar(x) > 45) x <- substr(x, 1, 45)

            x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],collapse=" "), "...", sep="")
            return(x)
        }
        else
        {
        return(x)
        }
}

labels=(sapply(
    levels(go_enrich_df$Description)[as.factor(go_enrich_df$Description)],
    shorten_names))

labels=(sapply( levels(factor(go_enrich_df$Description))[as.factor(go_enrich_df$number)], shorten_names))
names(labels) = rev(1:nrow(go_enrich_df))
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
library(ggplot2)

#go_enrich_df$number <- as.numeric(paste(1:40))
g <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
            geom_bar(stat="identity", width=0.8) + coord_flip() +
            scale_fill_manual(values = CPCOLS) + theme_bw() +
            xlab("GO term") +
            scale_x_discrete(labels=rev(go_enrich_df$Description)) +
            theme(axis.text=element_text(face = "bold", color="gray50")) +
            labs(title = "The Most Enriched GO Terms")

pdf('G1.GO_ALL.pdf',20,10)
g
dev.off()

g1 <- barplot(ego_BP,showCategory = 18,title="The GO_BP enrichment analysis of all DEGs ")+
              scale_size(range=c(2, 12))+
              scale_x_discrete(labels = function(ego_BP) str_wrap(ego_BP,width = 25))

pdf('G2.GO_BP.pdf')
g1
dev.off()

g2 <- barplot(ego_CC,showCategory = 18,title="The GO_CC enrichment analysis of all DEGs ")+
              scale_size(range=c(2, 12))+
              scale_x_discrete(labels = function(ego_CC) str_wrap(ego_CC,width = 25))

pdf('G3.GO_CC.pdf')
g2
dev.off()

g3 <- barplot(ego_MF,showCategory = 18,title="The GO_MF enrichment analysis of all DEGs ")+
              scale_size(range=c(2, 12))+
              scale_x_discrete(labels = function(ego_MF) str_wrap(ego_MF,width = 25))

pdf('G4.GO_BP.pdf')
g3
dev.off()

pdf('G5.MFplot.pdf')
plotGOgraph(ego_MF)
dev.off()

pdf('G6.CCplot.pdf')
plotGOgraph(ego_CC)
dev.off()

pdf('G7.BPplot.pdf')
plotGOgraph(ego_BP)
dev.off()


dir.create("/bigdata/wangzhang_guest/chenpeng_project/06_result/19_LJX_RNAseq/12.KEGG")
setwd("/bigdata/wangzhang_guest/chenpeng_project/06_result/19_LJX_RNAseq/12.KEGG")

kk<-enrichKEGG(gene =gene.df$ENTREZID,
               organism = 'mmu',
               pvalueCutoff = 0.05)

k1 <- barplot(kk,showCategory = 25, title="The KEGG enrichment analysis of all DEGs")
k1 + scale_size(range=c(2, 12))
k1 + scale_x_discrete(labels=function(kk) str_wrap(kk,width = 25))
k2 <- dotplot(kk,showCategory = 25,orderBy = "GeneRatio", title="The KEGG enrichment analysis of all DEGs")
k2 + scale_size(range=c(2, 12))
k2 + scale_x_discrete(labels=function(kk) str_wrap(kk,width = 25))

pdf('K1.KEGG_barplot.pdf')
k1
dev.off()

pdf('K2.KEGG_dotplot.pdf')
k2
dev.off()                      
