load("~/Pro_TNBC/paper/data/results/section_4/TNBC_BRCA1_0554_tumor_scRNA(tsne_by_cor).RData")
load("~/Pro_TNBC/paper/data/results/section_4/GSE202771_HDQP1_scRNA.RData")
GSE202771_HDQP1_scRNA                   <- HDQP1_scRNA
load("~/Pro_TNBC/paper/data/results/section_4/GSE173634_HDQP1_scRNA.RData")
GSE173634_HDQP1_scRNA                   <- HDQP1_scRNA

table(TNBC_BRCA1_0554_tumor_scRNA$subtype)
TNBC_BRCA1_0554_tumor_scRNA             <- subset(TNBC_BRCA1_0554_tumor_scRNA,subtype=="Basal"|subtype=="Claudin_low")
Idents(TNBC_BRCA1_0554_tumor_scRNA)     <- "subtype"
TNBC0554_pesudobulk                     <- AggregateExpression(
  TNBC_BRCA1_0554_tumor_scRNA,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "count",
  verbose = TRUE
)
TNBC0554_pesudobulk                     <- TNBC0554_pesudobulk$RNA
TNBC0554_pesudobulk_cpm                 <- apply(TNBC0554_pesudobulk ,2, function(x) { x/sum(x)*1000000 })
TNBC0554_pesudobulk_log2cpm             <- log2(TNBC0554_pesudobulk_cpm +1)
colnames(TNBC0554_pesudobulk_log2cpm)   <- c("TNBC0554_BL","TNBC0554_CL")


Idents(GSE173634_HDQP1_scRNA)           <- "subtype"
HDQP1_GSE173634_pesudobulk              <- AggregateExpression(
  GSE173634_HDQP1_scRNA,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "count",
  verbose = TRUE
)
HDQP1_GSE173634_pesudobulk                   <- HDQP1_GSE173634_pesudobulk$RNA
HDQP1_GSE173634_pesudobulk_cpm               <- apply(HDQP1_GSE173634_pesudobulk ,2, function(x) { x/sum(x)*1000000 })
HDQP1_GSE173634_pesudobulk_log2cpm           <- log2(HDQP1_GSE173634_pesudobulk_cpm +1)
colnames(HDQP1_GSE173634_pesudobulk_log2cpm) <- c("HDQP1_GSE173634_BL","HDQP1_GSE173634_CL")

table(GSE202771_HDQP1_scRNA$subtype)
GSE202771_HDQP1_scRNA                        <- subset(GSE202771_HDQP1_scRNA,subtype=="Basal"|subtype=="Claudin_low")
Idents(GSE202771_HDQP1_scRNA)                <- "subtype"
HDQP1_GSE202771_pesudobulk                   <- AggregateExpression(
  GSE202771_HDQP1_scRNA,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "count",
  verbose = TRUE
)
HDQP1_GSE202771_pesudobulk                   <- HDQP1_GSE202771_pesudobulk$RNA
HDQP1_GSE202771_pesudobulk_cpm               <- apply(HDQP1_GSE202771_pesudobulk ,2, function(x) { x/sum(x)*1000000 })
HDQP1_GSE202771_pesudobulk_log2cpm           <- log2(HDQP1_GSE202771_pesudobulk_cpm +1)
colnames(HDQP1_GSE202771_pesudobulk_log2cpm) <- c("HDQP1_GSE202771_BL","HDQP1_GSE202771_CL")
HDQP1_GSE202771_pesudobulk_log2cpm_df        <- as.data.frame(HDQP1_GSE202771_pesudobulk_log2cpm)
HDQP1_GSE202771_pesudobulk_log2cpm_df$SYMBOL <- rownames(HDQP1_GSE202771_pesudobulk_log2cpm)
symbol                                       <- bitr(HDQP1_GSE202771_pesudobulk_log2cpm_df$SYMBOL,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = "org.Hs.eg.db")
HDQP1_GSE202771_pesudobulk_log2cpm_df        <- merge(symbol,HDQP1_GSE202771_pesudobulk_log2cpm_df,by="SYMBOL")
HDQP1_GSE202771_pesudobulk_log2cpm_df        <- HDQP1_GSE202771_pesudobulk_log2cpm_df[!duplicated(HDQP1_GSE202771_pesudobulk_log2cpm_df$ENSEMBL),]
library(tibble)
rownames(HDQP1_GSE202771_pesudobulk_log2cpm_df) <- NULL
HDQP1_GSE202771_pesudobulk_log2cpm_df           <- column_to_rownames(HDQP1_GSE202771_pesudobulk_log2cpm_df,var = "ENSEMBL")
HDQP1_GSE202771_pesudobulk_log2cpm              <- as.matrix(HDQP1_GSE202771_pesudobulk_log2cpm_df[,-1])


####* Compute CCLE rna seq marker genes ##########
load("Pro_TNBC/data/CCLE/CCLE.RData")
CCLE.median                 <- apply(CCLE.log2.rpkm.matrix,1,median)
CCLE.expressed.gene         <- names(CCLE.median)[CCLE.median > 1]
tmp                         <- CCLE.log2.rpkm.matrix[CCLE.expressed.gene,]
tmp.rank                    <- apply(tmp,2,rank)
rank.mean                   <- apply(tmp.rank,1,mean)
rank.sd                     <- apply(tmp.rank,1,sd)
plot(x=rank.mean,y=rank.sd)
lowess(x=rank.mean,y=rank.sd) %>% lines(lwd=5,col='red')
CCLE.rna.seq.marker.gene.1000          <- names(sort(rank.sd,decreasing =TRUE))[1:1000]
save(CCLE.rna.seq.marker.gene.1000,file = "Pro_TNBC/output/data/CCLE/CCLE.rna.seq.marker.gene.1000.RData")

genes                                  <- list(rownames(TNBC0554_pesudobulk_log2cpm),rownames(HDQP1_GSE173634_pesudobulk_log2cpm),rownames(HDQP1_GSE202771_pesudobulk_log2cpm))
common_genes                           <- Reduce(intersect,genes)
TNBC0554_pesudobulk_log2cpm_CO         <- TNBC0554_pesudobulk_log2cpm[common_genes,]
HDQP1_GSE173634_pesudobulk_log2cpm_CO  <- HDQP1_GSE173634_pesudobulk_log2cpm[common_genes,]
HDQP1_GSE202771_pesudobulk_log2cpm_CO  <- HDQP1_GSE202771_pesudobulk_log2cpm[common_genes,]
identical(rownames(TNBC0554_pesudobulk_log2cpm_CO),rownames(HDQP1_GSE173634_pesudobulk_log2cpm_CO))
identical(rownames(TNBC0554_pesudobulk_log2cpm_CO),rownames(HDQP1_GSE202771_pesudobulk_log2cpm_CO))
log2cpm_co                             <- bind_cols(TNBC0554_pesudobulk_log2cpm_CO,HDQP1_GSE173634_pesudobulk_log2cpm_CO,HDQP1_GSE202771_pesudobulk_log2cpm_CO)
log2cpm_co                             <- as.matrix(log2cpm_co)
rownames(log2cpm_co)                   <- common_genes
BL_log2cpm_co                          <- log2cpm_co[,c(1,3,5)]
CL_log2cpm_co                          <- log2cpm_co[,c(2,4,6)]

####* Function to pick out cell line ##########
library(foreach)
pick.out.cell.line      <- function(expr.of.samples,expr.of.cell.lines,marker.gene){
  marker.gene           <- intersect(rownames(expr.of.samples),(marker.gene))
  marker.gene           <- intersect(rownames(expr.of.cell.lines),(marker.gene))
  correlation.matrix    <- cor(expr.of.samples[marker.gene,],expr.of.cell.lines[marker.gene,],method='spearman')
  cell.line.median.cor  <- apply(correlation.matrix,2,function(x) median(x[is.na(x) == FALSE])) %>% sort(decreasing = TRUE)
  best.cell.line        <- names(cell.line.median.cor)[1]
  p.value.vec           <- foreach(cell.line= setdiff(names(cell.line.median.cor),best.cell.line),.combine='c') %do% { #Find the different elements in vector x and vector y (only the different elements in x are taken)
    v                   <- correlation.matrix[,cell.line]
    p.value             <- wilcox.test(correlation.matrix[,best.cell.line],v,alternative = 'greater',paired = TRUE)$p.value
  }
  names(p.value.vec)    <- setdiff(names(cell.line.median.cor),best.cell.line)
  fdr.vec               <- p.adjust(p.value.vec,method='fdr')
  list(cell.line.median.cor=cell.line.median.cor,best.cell.line=best.cell.line,compare.fdr.vec=fdr.vec,correlation.matrix = correlation.matrix )
}


BL_rs                      <- pick.out.cell.line(expr.of.samples = BL_log2cpm_co,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
BL_median_cor              <- BL_rs$cell.line.median.cor %>% as.data.frame()
colnames(BL_median_cor)[1] <- "median_cor"
BL_median_cor$CCLE_Name    <- rownames(BL_median_cor)
BL_median_cor              <- merge(BL_median_cor,ccle_info[,c(4,13,21)],by="CCLE_Name")  

CL_rs                      <- pick.out.cell.line(expr.of.samples = CL_log2cpm_co,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
CL_median_cor              <- CL_rs$cell.line.median.cor %>% as.data.frame()
colnames(CL_median_cor)[1] <- "median_cor"
CL_median_cor$CCLE_Name    <- rownames(CL_median_cor)
CL_median_cor              <- merge(CL_median_cor,ccle_info[,c(4,13,21)],by="CCLE_Name")  

save(CL_rs,BL_rs,CL_median_cor,BL_median_cor,file = "Pro_TNBC/paper/data/results/section_4/cor.with.CCLE.RData")
