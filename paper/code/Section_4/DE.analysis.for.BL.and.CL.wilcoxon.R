
################################################################################
####Find the differentiated genes between Basal and Claudin-low cells
################################################################################

####*DE analysis of basal and claudin-low cells in TNBC 0554####
load("~/Pro_TNBC/paper/data/results/section_4/TNBC_0554/TNBC_BRCA1_0554_tumor_scRNA(tsne_by_cor).RData")
library(Seurat)
library(clusterProfiler)
Idents(TNBC_BRCA1_0554_tumor_scRNA)         <- "subtype"
de_results_TNBC0554_wilcox                  <- FindMarkers(TNBC_BRCA1_0554_tumor_scRNA, slot="data",ident.1 = "Claudin_low", ident.2 = "Basal",test.use = "wilcox",logfc.threshold = 0,pseudocount.use = 10^(-9))
de_results_TNBC0554_wilcox                  <- within(de_results_TNBC0554_wilcox,{
  sig                                       <- NA
  sig[p_val_adj < 0.05 &avg_log2FC > 0.5]   <- "up"
  sig[p_val_adj < 0.05 &avg_log2FC < -0.5]  <- "down"
})
de_results_TNBC0554_wilcox[is.na(de_results_TNBC0554_wilcox$sig),6]      <- "none"
table(de_results_TNBC0554_wilcox$sig)
de_results_TNBC0554_wilcox$ENSEMBL               <- rownames(de_results_TNBC0554_wilcox)
symbol                                           <- bitr(de_results_TNBC0554_wilcox$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
de_results_TNBC0554_wilcox                       <- merge(symbol,de_results_TNBC0554_wilcox,by="ENSEMBL")


####*HDQP1 in GSE173634####
####*
Idents(GSE173634.HDQP1.scRNA)                   <- "subtype"
de_results_HDQP1_GSE173634_wilcox               <- FindMarkers(GSE173634.HDQP1.scRNA,slot="data",ident.1 = "Claudin_low",ident.2 = "Basal",test.use = "wilcox",logfc.threshold = 0,pseudocount.use = 10^(-9))
de_results_HDQP1_GSE173634_wilcox               <- within(de_results_HDQP1_GSE173634_wilcox,{
  sig                                           <- NA
  sig[p_val_adj < 0.05 &avg_log2FC >0.5]        <- "up"
  sig[p_val_adj < 0.05 &avg_log2FC < -0.5]      <- "down"
})
de_results_HDQP1_GSE173634_wilcox[is.na(de_results_HDQP1_GSE173634_wilcox$sig),6] <- "none"
table(de_results_HDQP1_GSE173634_wilcox$sig)
de_results_HDQP1_GSE173634_wilcox$ENSEMBL               <- rownames(de_results_HDQP1_GSE173634_wilcox)
symbol                                                  <- bitr(de_results_HDQP1_GSE173634_wilcox$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
de_results_HDQP1_GSE173634_wilcox                       <- merge(symbol,de_results_HDQP1_GSE173634_wilcox,by="ENSEMBL")



####*HDQP1 in GSE202271####
load("~/Pro_TNBC/paper/data/results/section_4/GSE202771_HDQP1/GSE202771_HDQP1_scRNA.RData")
Idents(GSE202771_HDQP1_scRNA)                   <- "subtype"
de_results_HDQP1_GSE202771_wilcox               <- FindMarkers(GSE202771_HDQP1_scRNA,slot="data",ident.1 = "Claudin_low",ident.2 = "Basal",test.use = "wilcox",logfc.threshold = 0,pseudocount.use = 10^(-9))
de_results_HDQP1_GSE202771_wilcox               <- within(de_results_HDQP1_GSE202771_wilcox,{
  sig                                           <- NA
  sig[p_val_adj < 0.05 &avg_log2FC >0.5]        <- "up"
  sig[p_val_adj < 0.05 &avg_log2FC < -0.5]      <- "down"
})
de_results_HDQP1_GSE202771_wilcox[is.na(de_results_HDQP1_GSE202771_wilcox$sig),6] <- "none"
table(de_results_HDQP1_GSE202771_wilcox$sig)
de_results_HDQP1_GSE202771_wilcox$SYMBOL        <- rownames(de_results_HDQP1_GSE202771_wilcox)
library(clusterProfiler)
symbol                                                  <- bitr(de_results_HDQP1_GSE202771_wilcox$SYMBOL,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = "org.Hs.eg.db")
de_results_HDQP1_GSE202771_wilcox                       <- merge(de_results_HDQP1_GSE202771_wilcox,symbol,by="SYMBOL")


####*intersected de genes####
de_HDQP1_claudinlow_basal_intersect                 <- merge(de_results_HDQP1_GSE173634_wilcox,de_results_HDQP1_GSE202771_wilcox,by="ENSEMBL")
de_HDQP1_claudinlow_basal_intersect                 <- de_HDQP1_claudinlow_basal_intersect[,c(1,2,4,7,8,11,14,15)]
colnames(de_HDQP1_claudinlow_basal_intersect)[2:8]  <- c("SYMBOL","log2FC_HDQP1_GSE173634","padj_HDQP1_GSE173634","sig_1","log2FC_HDQP1_GSE202771","padj_HDQP1_GSE202771","sig_2")
de_intersect_gene                     <- merge(de_HDQP1_claudinlow_basal_intersect,de_results_TNBC0554_wilcox,by="ENSEMBL")
de_intersect_gene                     <- de_intersect_gene[,c(1:8,11,14,15)]
colnames(de_intersect_gene)[c(9:11)]  <- c("log2FC_TNBC0554","padj_TNBC0554","sig_0554")
de_intersect_gene                     <- within(de_intersect_gene,{
  sig                                                <- NA
  sig[sig_1 =="up"&sig_2=="up"&sig_0554=="up"]       <- "up"
  sig[sig_1 =="down"&sig_2=="down"&sig_0554=="down"] <- "down"
})
table(de_intersect_gene$sig)
save(de_intersect_gene,file = "Pro_TNBC/paper/data/results/section_5/de_intersect_gene_wilcox.RData")
save(de_results_TNBC0554_wilcox,file = "Pro_TNBC/paper/data/results/section_4/TNBC_0554/de_results_TNBC0554_wilcox.RData")
save(de_results_HDQP1_GSE173634_wilcox,file = "Pro_TNBC/paper/data/results/section_4/GSE173634_HDQP1/de_results_HDQP1_GSE173634_wilcox.RData")
save(de_results_HDQP1_GSE202771_wilcox,file = "Pro_TNBC/paper/data/results/section_4/GSE202771_HDQP1/de_results_HDQP1_GSE202771_wilcox.RData")

de_results_HDQP1_GSE173634_wilcox   <- subset(de_results_HDQP1_GSE173634_wilcox,sig=="up"|sig=="down")
de_results_HDQP1_GSE202771_wilcox   <- subset(de_results_HDQP1_GSE202771_wilcox,sig=="up"|sig=="down")
de_results_TNBC0554_wilcox          <- subset(de_results_TNBC0554_wilcox,sig=="up"|sig=="down")
de_intersect_gene_wilcox            <- de_intersect_gene[!is.na(de_intersect_gene$sig),]
de_results_HDQP1_GSE202771_wilcox   <- de_results_HDQP1_GSE202771_wilcox[!duplicated(de_results_HDQP1_GSE202771_wilcox$SYMBOL),]
write.csv(de_results_HDQP1_GSE202771_wilcox,file = "Pro_TNBC/paper/data/results/section_4/GSE202771_HDQP1/de_results_HDQP1_GSE202771_wilcox.csv")
write.csv(de_results_HDQP1_GSE173634_wilcox,file = "Pro_TNBC/paper/data/results/section_4/GSE173634_HDQP1/de_results_HDQP1_GSE173634_wilcox.csv")
write.csv(de_results_TNBC0554_wilcox,file = "Pro_TNBC/paper/data/results/section_4/TNBC_0554/de_results_TNBC0554_wilcox.csv")
write.csv(de_intersect_gene_wilcox,file = "Pro_TNBC/paper/data/results/section_5/de_intersect_gene.wilcox.csv")
