

####run seurat pipeline for HCC70_SCR####
library(Seurat)
dir=("/home/BioklabData/lijing_Data/KDELF3/CellRanger_V7/HCC70_SCR/outs/filtered_feature_bc_matrix")
HCC70_SCR_count                     <- Read10X(
  data.dir=dir,
  gene.column = 1,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE)
HCC70_SCR_scRNA                     <- CreateSeuratObject(counts = HCC70_SCR_count,min.cells = 3, min.features = 200)
#calculate the proportion of ribosomal genes in cells
HCC70_SCR_scRNA[["percent.mt"]]     <- PercentageFeatureSet(HCC70_SCR_scRNA,pattern = "^MT-")
head(HCC70_SCR_scRNA@meta.data,5)
VlnPlot(HCC70_SCR_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
HCC70_SCR_scRNA                     <- subset(HCC70_SCR_scRNA,nFeature_RNA >3000)

#normalized
HCC70_SCR_scRNA                     <- NormalizeData(HCC70_SCR_scRNA)

#features choosing
HCC70_SCR_scRNA                     <- FindVariableFeatures(HCC70_SCR_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                               <- head(VariableFeatures(HCC70_SCR_scRNA),10)

#data scaling 
all_gene                            <- rownames(HCC70_SCR_scRNA)
HCC70_SCR_scRNA                     <- ScaleData(HCC70_SCR_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data

#linear dimensionality reduction
HCC70_SCR_scRNA                     <- RunPCA(HCC70_SCR_scRNA,features = VariableFeatures(object=HCC70_SCR_scRNA))

#dimension selection
HCC70_SCR_scRNA                     <- JackStraw(HCC70_SCR_scRNA,num.replicate = 100)
HCC70_SCR_scRNA                     <- ScoreJackStraw(HCC70_SCR_scRNA,dims = 1:20)
JackStrawPlot(HCC70_SCR_scRNA,dims=1:20)
ElbowPlot(HCC70_SCR_scRNA,ndims = 50)#choose dims of 40
#cell clustering
HCC70_SCR_scRNA                     <- FindNeighbors(HCC70_SCR_scRNA,dims = 1:40)
HCC70_SCR_scRNA                     <- FindClusters(HCC70_SCR_scRNA,resolution = 0.1)
table(Idents(HCC70_SCR_scRNA))

#unlinear dimension reduction(UMAP or tSNE)
HCC70_SCR_scRNA                     <- RunUMAP(HCC70_SCR_scRNA,dims = 1:40)
HCC70_SCR_scRNA                     <- RunTSNE(HCC70_SCR_scRNA,dims = 1:40,check_duplicates = FALSE)
DimPlot(HCC70_SCR_scRNA,reduction = "umap",label = T)
DimPlot(HCC70_SCR_scRNA,reduction = "tsne",label = T)

####UBS93 subtyping for HCC70_SCR####
library(tibble)
load("Pro_TNBC/paper/data/method/UBS93.data.RData")
HCC70_SCR_mat               <- HCC70_SCR_scRNA@assays$RNA@data %>% as.matrix()
predictor.res               <- breast.cancer.predictor(expr.of.sample = HCC70_SCR_mat,
                                                       expr.of.centroid = UBS93.data$UBS93.centroid,
                                                       marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                       HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                       ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
HCC70_SCR_subtype           <- predictor.res$subtype
table(HCC70_SCR_subtype)
HCC70_SCR_metadata          <- HCC70_SCR_scRNA@meta.data
identical(rownames(HCC70_SCR_metadata),rownames(HCC70_SCR_subtype))
HCC70_SCR_metadata$subtype  <- HCC70_SCR_subtype$subtype
HCC70_SCR_scRNA@meta.data   <- HCC70_SCR_metadata
Idents(HCC70_SCR_scRNA)     <- "subtype"
DimPlot(HCC70_SCR_scRNA,reduction = "umap",label = T)
DimPlot(HCC70_SCR_scRNA,reduction = "tsne",label = T)
saveRDS(HCC70_SCR_scRNA,file = "Pro_TNBC/paper/data/results/section_5/HCC70_SCR_scRNA.rds")
HCC70_SCR_cor               <- predictor.res$cor.matrix %>% as.data.frame()
HCC70_SCR_cor$subtype       <- predictor.res$subtype$subtype


HCC70_SCR_exp_df            <- as.data.frame(HCC70_SCR_mat)
HCC70_SCR_exp_df$ENSEMBL    <- rownames(HCC70_SCR_exp_df) 
symbol                      <- bitr(HCC70_SCR_exp_df$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb="org.Hs.eg.db")
HCC70_SCR_exp_df            <- merge(symbol,HCC70_SCR_exp_df,by="ENSEMBL")
marker.gene                 <- c("ELF3","CLDN3","CLDN4","CLDN7","EPCAM","VIM")
HCC70_SCR_exp_df            <- HCC70_SCR_exp_df[HCC70_SCR_exp_df$SYMBOL %in% marker.gene,]
rownames(HCC70_SCR_exp_df)  <- NULL
HCC70_SCR_exp_df            <- column_to_rownames(HCC70_SCR_exp_df,var = "SYMBOL")
HCC70_SCR_exp_df            <- HCC70_SCR_exp_df[,-1]
HCC70_SCR_exp_df            <- as.matrix(HCC70_SCR_exp_df) %>% t() %>% as.data.frame()
HCC70_SCR_exp_df$subtype    <- HCC70_SCR_scRNA$subtype
ggplot(HCC70_SCR_exp_df,aes(x=subtype,y=VIM))+
  geom_boxplot()+
  geom_point(size=1,position="jitter") 
kruskal.test(ELF3 ~subtype,data = HCC70_SCR_exp_df)#p-value = 1.83e-14
saveRDS(HCC70_SCR_exp_df,"Pro_TNBC/paper/data/results/section_5/HCC70_SCR_exp_df.rds")

####run seurat pipeline for HCC70_SHELF3####
library(Seurat)
dir=("/home/BioklabData/lijing_Data/KDELF3/CellRanger_V7/HCC70_SHELF3/outs/filtered_feature_bc_matrix/")
HCC70_SHELF3_count                     <- Read10X(
  data.dir=dir,
  gene.column = 1,
  cell.column = 1,
  unique.features = TRUE,
  strip.suffix = FALSE)
HCC70_SHELF3_scRNA                    <- CreateSeuratObject(counts = HCC70_SHELF3_count,min.cells = 3, min.features = 200)
#calculate the proportion of ribosomal genes in cells
HCC70_SHELF3_scRNA[["percent.mt"]]     <- PercentageFeatureSet(HCC70_SHELF3_scRNA,pattern = "^MT-")
head(HCC70_SHELF3_scRNA@meta.data,5)
VlnPlot(HCC70_SHELF3_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
HCC70_SHELF3_scRNA                     <- subset(HCC70_SHELF3_scRNA,nFeature_RNA >0 )

#normalized
HCC70_SHELF3_scRNA                     <- NormalizeData(HCC70_SHELF3_scRNA)

#features choosing
HCC70_SHELF3_scRNA                     <- FindVariableFeatures(HCC70_SHELF3_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                               <- head(VariableFeatures(HCC70_SHELF3_scRNA),10)

#data scaling 
all_gene                               <- rownames(HCC70_SHELF3_scRNA)
HCC70_SHELF3_scRNA                     <- ScaleData(HCC70_SHELF3_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data

#linear dimensionality reduction
HCC70_SHELF3_scRNA                     <- RunPCA(HCC70_SHELF3_scRNA,features = VariableFeatures(object=HCC70_SHELF3_scRNA))

#dimension selection
HCC70_SHELF3_scRNA                     <- JackStraw(HCC70_SHELF3_scRNA,num.replicate = 100)
HCC70_SHELF3_scRNA                     <- ScoreJackStraw(HCC70_SHELF3_scRNA,dims = 1:20)
JackStrawPlot(HCC70_SHELF3_scRNA,dims=1:20)
ElbowPlot(HCC70_SHELF3_scRNA,ndims = 50)#choose dims of 30
#cell clustering
HCC70_SHELF3_scRNA                     <- FindNeighbors(HCC70_SHELF3_scRNA,dims = 1:40)
HCC70_SHELF3_scRNA                     <- FindClusters(HCC70_SHELF3_scRNA,resolution = 0.1)
table(Idents(HCC70_SHELF3_scRNA))

#unlinear dimension reduction(UMAP or tSNE)
HCC70_SHELF3_scRNA                     <- RunUMAP(HCC70_SHELF3_scRNA,dims = 1:40)
HCC70_SHELF3_scRNA                     <- RunTSNE(HCC70_SHELF3_scRNA,dims = 1:40,check_duplicates = FALSE)
DimPlot(HCC70_SHELF3_scRNA,reduction = "umap",label = T)
DimPlot(HCC70_SHELF3_scRNA,reduction = "tsne",label = T)

####UBS93 subtyping for HCC70_SHELF3####
library(tibble)
load("Pro_TNBC/paper/data/method/UBS93.data.RData")
HCC70_SHELF3_mat               <- HCC70_SHELF3_scRNA@assays$RNA@data %>% as.matrix()
predictor.res                  <- breast.cancer.predictor(expr.of.sample = HCC70_SHELF3_mat,
                                                          expr.of.centroid = UBS93.data$UBS93.centroid,
                                                          marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                          HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                          ssgsea.cutoff = UBS93.data$ssgsea.cutoff)
saveRDS(predictor.res,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_SHELF3_UBS93_predictionres.rds")
HCC70_SHELF3_subtype           <- predictor.res$subtype
table(HCC70_SHELF3_subtype)
HCC70_SHELF3_metadata          <- HCC70_SHELF3_scRNA@meta.data
identical(rownames(HCC70_SHELF3_metadata),rownames(HCC70_SHELF3_subtype))
HCC70_SHELF3_metadata$subtype  <- HCC70_SHELF3_subtype$subtype
HCC70_SHELF3_scRNA@meta.data   <- HCC70_SHELF3_metadata
Idents(HCC70_SHELF3_scRNA)     <- "subtype"
DimPlot(HCC70_SHELF3_scRNA,reduction = "umap",label = T)
DimPlot(HCC70_SHELF3_scRNA,reduction = "tsne",label = T)
saveRDS(HCC70_SHELF3_scRNA,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_SHELF3_scRNA.rds")
HCC70_SHELF3_cor               <- predictor.res$cor.matrix %>% as.data.frame()
HCC70_SHELF3_cor$subtype       <- predictor.res$subtype$subtype

####*TSNE####
HCC70_SHELF3_exp                       <- as.matrix(HCC70_SHELF3_scRNA@assays$RNA@data)
HCC70_SHELF3_exp_UBS83                 <- HCC70_SHELF3_exp[rownames(HCC70_SHELF3_exp) %in% UBS93.data$UBS93.gene.df$ENSEMBL,]
library(Rtsne)
cor_matrix                                      <- cor(HCC70_SHELF3_exp_UBS83,method = "spearman")
dissimilarity_matrix                            <- 1 - cor_matrix
normalized_dissimilarity_matrix                 <- dissimilarity_matrix / max(dissimilarity_matrix)
HCC70_SHELF3_scRNA                     <- RunTSNE(
  HCC70_SHELF3_scRNA,
  reduction = "pca",
  cells = NULL,
  dims = 1:40,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = normalized_dissimilarity_matrix,
  reduction.name = "tsne",
  reduction.key = "tSNE_"
)


####* marker genes in tsne plot####
p1 <- FeaturePlot(HCC70_SHELF3_scRNA,features = c("ENSG00000163435"),#ELF3
                  reduction = "tsne",pt.size = 3)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5))+ggtitle("ELF3")+theme_bw(base_size = 55)+
  xlab("")+ylab("") 

p2 <- FeaturePlot(HCC70_SHELF3_scRNA,features = c("ENSG00000165215"),#CLDN3
                  reduction = "tsne",pt.size = 3)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5))+ggtitle("CLDN3")+theme_bw(base_size = 55)+
  xlab("")+ylab("") 
p3 <- FeaturePlot(HCC70_SHELF3_scRNA,features = c("ENSG00000189143"),#CLDN4
                  reduction = "tsne",pt.size = 3)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5))+ggtitle("CLDN4")+theme_bw(base_size = 55)+
  xlab("")+ylab("")
p4 <- FeaturePlot(HCC70_SHELF3_scRNA,features = c("ENSG00000181885"),#CLDN7
                  reduction = "tsne",pt.size = 3)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5))+ggtitle("CLDN7")+theme_bw(base_size = 55) +
  xlab("")+ylab("")
p6 <- FeaturePlot(HCC70_SHELF3_scRNA,features = c("ENSG00000026025"),#VIM
                  reduction = "tsne",pt.size = 3)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5))+ggtitle("VIM")+theme_bw(base_size = 55) +
  xlab("")+ylab("")
p5 <- FeaturePlot(HCC70_SHELF3_scRNA,features = c("ENSG00000119888"),#EPCAM
                  reduction = "tsne",pt.size = 3)+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5))+ggtitle("EPCAM")+theme_bw(base_size = 55) +
  xlab("")+ylab("")
Idents(HCC70_SHELF3_scRNA)   <- "subtype"
cols                <- c("grey","blue","red","yellow")
p7                  <- DimPlot(HCC70_SHELF3_scRNA,reduction = "tsne",label = F,pt.size = 3,label.size =20,cols=cols )+
  theme( 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.key.size = unit(25, "pt"),
    legend.key.height = unit(25, "pt"),
    legend.key.width = unit(35, "pt"),plot.margin=unit(rep(1,4),'cm')
  )+theme_bw(base_size = 55)
combined_plot       <- plot_grid(p1, p2,p3,p4,p5,p6 ,p7,ncol = 3) +theme_bw(base_size = 55)
ggsave(combined_plot,filename = "Pro_TNBC/paper/plot/section_5/TSNE.HCC70.SHELF3.new.pdf",width = 50,height = 50,limitsize = F)

saveRDS(HCC70_SHELF3_scRNA,file = "Pro_TNBC/paper/data/results/section_5/HCC70_SHELF3_scRNA.rds")

####* correlation coefficients in tsne plot####
HCC70_SHELF3_UBS93_predictionres <- readRDS("~/Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_SHELF3_UBS93_predictionres.rds")
HCC70_SHELF3_scRNA@meta.data     <- cbind(HCC70_SHELF3_scRNA@meta.data,HCC70_SHELF3_UBS93_predictionres$cor.matrix)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
p1 <- FeaturePlot(HCC70_SHELF3_scRNA,features = c("Basal"),#Basal
                  reduction = "tsne",pt.size = 3)+scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5))+ggtitle("Basal")+theme_bw(base_size = 55) +
  xlab("")+ylab("")
p2 <- FeaturePlot(HCC70_SHELF3_scRNA,features = c("Claudin_low"),#EPCAM
                  reduction = "tsne",pt.size = 3)+scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5))+ggtitle("Claudin_low")+theme_bw(base_size = 55) +
  xlab("")+ylab("")
p3 <- FeaturePlot(HCC70_SHELF3_scRNA,features = c("H_or_L"),#EPCAM
                  reduction = "tsne",pt.size = 3)+scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))+
  theme( 
    panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
    axis.ticks = element_blank(),axis.text.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "none", 
    plot.title = element_text(size=55,hjust=0.5))+ggtitle("H_or_L")+theme_bw(base_size = 55) +
  xlab("")+ylab("")
combined_plot       <- plot_grid(p1, p2,p3,ncol = 3) +theme_bw(base_size = 55)
ggsave(combined_plot,filename = "Pro_TNBC/paper/plot/section_5/TSNE.HCC70.SHELF3.cor.new.pdf",width = 50,height = 15,limitsize = F)

####boxplot####
HCC70_SHELF3_exp_df            <- as.data.frame(HCC70_SHELF3_mat)
HCC70_SHELF3_exp_df$ENSEMBL    <- rownames(HCC70_SHELF3_exp_df) 
symbol                         <- bitr(HCC70_SHELF3_exp_df$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb="org.Hs.eg.db")
HCC70_SHELF3_exp_df            <- merge(symbol,HCC70_SHELF3_exp_df,by="ENSEMBL")
marker.gene                    <- c("ELF3","CLDN3","CLDN4","CLDN7","EPCAM","VIM")
HCC70_SHELF3_exp_df            <- HCC70_SHELF3_exp_df[HCC70_SHELF3_exp_df$SYMBOL %in% marker.gene,]
rownames(HCC70_SHELF3_exp_df)  <- NULL
HCC70_SHELF3_exp_df            <- column_to_rownames(HCC70_SHELF3_exp_df,var = "SYMBOL")
HCC70_SHELF3_exp_df            <- HCC70_SHELF3_exp_df[,-1]
HCC70_SHELF3_exp_df            <- as.matrix(HCC70_SHELF3_exp_df) %>% t() %>% as.data.frame()
HCC70_SHELF3_exp_df$subtype    <- HCC70_SHELF3_scRNA$subtype
ggplot(HCC70_SHELF3_exp_df,aes(x=subtype,y=ELF3))+
  geom_boxplot()+
  geom_point(size=1,position="jitter") 
kruskal.test(ELF3 ~subtype,data = HCC70_SHELF3_exp_df)#p-value = 1.83e-14
saveRDS(HCC70_SHELF3_exp_df,"Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_SHELF3_exp_df.rds")

p1                  <- ggplot(HCC70_SHELF3_exp_df,aes(x=subtype,y=ELF3))+
  geom_boxplot() +stat_compare_means(method = "kruskal.test")
p2                  <- ggplot(HCC70_SHELF3_exp_df,aes(x=subtype,y=CLDN3))+
  geom_boxplot() +stat_compare_means(method = "kruskal.test")
p3                  <- ggplot(HCC70_SHELF3_exp_df,aes(x=subtype,y=CLDN4))+
  geom_boxplot() +stat_compare_means(method = "kruskal.test")
p4                  <- ggplot(HCC70_SHELF3_exp_df,aes(x=subtype,y=CLDN7))+
  geom_boxplot() +stat_compare_means(method = "kruskal.test")
p5                  <- ggplot(HCC70_SHELF3_exp_df,aes(x=subtype,y=EPCAM))+
  geom_boxplot() +stat_compare_means(method = "kruskal.test")
p6                  <- ggplot(HCC70_SHELF3_exp_df,aes(x=subtype,y=VIM))+
  geom_boxplot() +stat_compare_means(method = "kruskal.test")
combined_plot       <- plot_grid(p1, p2,p3,p4,p5,p6, ncol = 3) 
combined_plot


library(ggplot2)
library(ggpubr)
library(cowplot)
HCC70_exp_df        <- rbind(HCC70_SCR_exp_df,HCC70_SHELF3_exp_df)
HCC70_exp_df$group  <- c(rep("HCC70_SCR",6161),rep("HCC70_SHELF3",10386))
p1                  <- ggplot(HCC70_exp_df,aes(x=group,y=ELF3))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")
p2                  <- ggplot(HCC70_exp_df,aes(x=group,y=CLDN3))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")
p3                  <- ggplot(HCC70_exp_df,aes(x=group,y=CLDN4))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")
p4                  <- ggplot(HCC70_exp_df,aes(x=group,y=CLDN7))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")
p5                  <- ggplot(HCC70_exp_df,aes(x=group,y=EPCAM))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")
p6                  <- ggplot(HCC70_exp_df,aes(x=group,y=VIM))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")
combined_plot       <- plot_grid(p1, p2,p3,p4,p5,p6, ncol = 3) 
combined_plot


HCC70_cor        <- rbind(HCC70_SCR_cor,HCC70_SHELF3_cor)
HCC70_cor$group  <- c(rep("HCC70_SCR",6161),rep("HCC70_SHELF3",10386))
HCC70_cor_basal  <- subset(HCC70_cor,subtype=="Basal")
ggplot(HCC70_cor_basal,aes(x=group,y=Basal))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")
ggplot(HCC70_cor,aes(x=group,y=Basal))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")
ggplot(HCC70_cor,aes(x=group,y=Claudin_low))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")


####Pseudo bulk by sum####
Idents(HCC70_SCR_scRNA)     <- "orig.ident"
Idents(HCC70_SHELF3_scRNA)  <- "orig.ident"
HCC70_SCR_pesudobulk        <- AggregateExpression(
  HCC70_SCR_scRNA,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "count",
  verbose = TRUE
)
HCC70_SHELF3_pesudobulk      <- AggregateExpression(
  HCC70_SHELF3_scRNA,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "count",
  verbose = TRUE
)
HCC70_SCR_pesudobulk         <- HCC70_SCR_pesudobulk$RNA
library(edgeR)
HCC70_SCR_pesudobulk_cpm     <- cpm(HCC70_SCR_pesudobulk)
HCC70_SCR_pesudobulk_cpm["ENSG00000163435",]#ELF3
HCC70_SCR_pesudobulk_cpm["ENSG00000165215",]#CLDN3

HCC70_SHELF3_pesudobulk      <- HCC70_SHELF3_pesudobulk$RNA
library(edgeR)
HCC70_SHELF3_pesudobulk_cpm  <- cpm(HCC70_SHELF3_pesudobulk)
HCC70_SHELF3_pesudobulk_cpm["ENSG00000163435",]
HCC70_SHELF3_pesudobulk_cpm["ENSG00000165215",]

HCC70_SCR_pesudobulk_sum       <- apply(HCC70_SCR_count,1,sum)
HCC70_SCR_pesudobulk_sum_cpm   <- cpm(HCC70_SCR_pesudobulk_sum)
HCC70_SCR_pesudobulk_sum_cpm["ENSG00000163435",]#ELF3

#average
HCC70_SCR_count                    <- as.matrix(HCC70_SCR_scRNA@assays$RNA@counts) 
HCC70_SCR_pesudobulk_average       <- apply(HCC70_SCR_count,1,mean)
HCC70_SCR_pesudobulk_average_cpm   <- cpm(HCC70_SCR_pesudobulk_average)
HCC70_SCR_pesudobulk_average_cpm["ENSG00000163435",]#ELF3


HCC70_SHELF3_count                    <- as.matrix(HCC70_SHELF3_scRNA@assays$RNA@counts) 
HCC70_SHELF3_pesudobulk_average       <- apply(HCC70_SHELF3_count,1,mean)
HCC70_SHELF3_pesudobulk_average_cpm   <- cpm(HCC70_SHELF3_pesudobulk_average)
HCC70_SHELF3_pesudobulk_average_cpm["ENSG00000163435",]#ELF3

####DE analysis for ELF3 high and low in SHELF3####
HCC70_SHELF3_cpm                 <- HCC70_SHELF3_scRNA@assays$RNA@data %>% as.matrix()%>% t() %>% as.data.frame()
HCC70_SHELF3_cpm                 <- HCC70_SHELF3_cpm[order(HCC70_SHELF3_cpm$ENSG00000163435,decreasing = F),]
all(diff(HCC70_SHELF3_cpm[,"ENSG00000163435"]) >= 0)
HCC70_SHELF3_cpm_high_low        <- HCC70_SHELF3_cpm[c(1:500,9887:10386 ),]
HCC70_SHELF3_high_low_scRNA      <- HCC70_SHELF3_scRNA[,rownames(HCC70_SHELF3_cpm_high_low)]
HCC70_SHELF3_high_low_scRNA$group    <- c(rep("ELF3_low",500),rep("ELF3_high",500))


Idents(HCC70_SHELF3_high_low_scRNA)  <- "group"
#wilcox
HCC70_SHELF3_high_low_de_res            <- FindMarkers(HCC70_SHELF3_high_low_scRNA,slot = "data",test.use = "wilcox",ident.1="ELF3_low",ident.2="ELF3_high",logfc.threshold = 0)
HCC70_SHELF3_high_low_deres             <- within(HCC70_SHELF3_high_low_de_res,{
  sig                                   <- NA
  sig[p_val_adj < 0.05 & avg_log2FC > 0]   <- "up"
  sig[p_val_adj < 0.05 & avg_log2FC < 0]  <- "down"
})
HCC70_SHELF3_high_low_deres[is.na(HCC70_SHELF3_high_low_deres$sig),6] <- "none"
#
symbol                                <- bitr(rownames(HCC70_SHELF3_high_low_deres),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
HCC70_SHELF3_high_low_deres$ENSEMBL   <- rownames(HCC70_SHELF3_high_low_deres)
HCC70_SHELF3_high_low_deres_symbol    <- merge(HCC70_SHELF3_high_low_deres,symbol,by="ENSEMBL")
saveRDS(HCC70_SHELF3_high_low_deres_symbol,file = "Pro_TNBC/paper/data/results/section_5/HCC70_SHELF3_high_low_deres_symbol.rds")
HCC70_SHELF3_low_VS_high_deres_symbol <- HCC70_SHELF3_high_low_deres_symbol
write.csv(HCC70_SHELF3_low_VS_high_deres_symbol,file = "Pro_TNBC/paper/data/results/section_5/HCC70_SHELF3_low_VS_high_deres_symbol.csv")
saveRDS(HCC70_SHELF3_high_low_scRNA,file = "Pro_TNBC/paper/data/results/section_5/HCC70_SHELF3_high_low_scRNA.rds")

HCC70_SHELF3_high_low_exp_df        <- HCC70_SHELF3_high_low_scRNA@assays$RNA@data %>% as.matrix() %>% t()%>% as.data.frame()
HCC70_SHELF3_high_low_exp_df$group  <- c(rep("ELF3_low",500),rep("ELF3_high",500))
p1   <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000165215))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("CLDN3")
p2  <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000189143))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("CLDN4")
p3 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000181885))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("CLDN7")
p4 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000119888))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("EPCAM")
p5 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000026025))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("VIM")
p6 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000163435))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("ELF3")
p7 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000186847))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("KRT14")
p8 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000186081))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("KRT5")
combined_plot       <- plot_grid(p1, p2,p3 ,p4,p5,p6,p7,p8,ncol = 3) 
combined_plot


####correlation coefficients for ELF3 high and low in SHELF3####
HCC70_SHELF3_exp_df_high_low_cor          <- HCC70_SHELF3_cor[rownames(HCC70_SHELF3_high_low_exp_df),]
identical(rownames(HCC70_SHELF3_exp_df_high_low_cor),rownames(HCC70_SHELF3_high_low_exp_df))
HCC70_SHELF3_exp_df_high_low_cor$group    <- c(rep("ELF3_low",500),rep("ELF3_high",500))
p1                  <- ggplot(HCC70_SHELF3_exp_df_high_low_cor,aes(x=group,y=Basal))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")
p2                  <- ggplot(HCC70_SHELF3_exp_df_high_low_cor,aes(x=group,y=Claudin_low))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")
p3 <-   ggplot(HCC70_SHELF3_exp_df_high_low_cor,aes(x=group,y=H_or_L))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")
combined_plot       <- plot_grid(p1, p2,p3 ,ncol = 3) 
combined_plot



####get pseudo bulk by the rank of ELF3 expression####
HCC70_SHELF3_cpm                 <- HCC70_SHELF3_scRNA@assays$RNA@data %>% as.matrix()%>% t() %>% as.data.frame()
HCC70_SHELF3_cpm                 <- HCC70_SHELF3_cpm[order(HCC70_SHELF3_cpm$ENSG00000163435,decreasing = F),]
all(diff(HCC70_SHELF3_cpm[,"ENSG00000163435"]) >= 0)
HCC70_SHELF3_count               <- HCC70_SHELF3_scRNA@assays$RNA@counts %>% as.matrix() %>% t() %>% as.data.frame()
HCC70_SHELF3_count               <- HCC70_SHELF3_count[rownames(HCC70_SHELF3_cpm),]
HCC70_SHELF3_count_sum                    <- apply(HCC70_SHELF3_count, 2, function(x) {
  matrix(x, nrow = 6, byrow  = F) %>% colSums()
})
rownames(HCC70_SHELF3_count_sum)          <- paste0("sample",1:1731)
library(edgeR)
HCC70_SHELF3_count_sum_log2cpm            <- apply(HCC70_SHELF3_count_sum,1,function(x){cpm(x,log=T)}) 
rownames(HCC70_SHELF3_count_sum_log2cpm)  <- colnames(HCC70_SHELF3_count_sum)

saveRDS(HCC70_SHELF3_count_sum_log2cpm,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_SHELF3_pusdobulk_log2cpm.rds")
load("Pro_TNBC/paper/data/method/UBS93.data.RData")
predictor.res                  <- breast.cancer.predictor(expr.of.sample = HCC70_SHELF3_count_sum_log2cpm,
                                                          expr.of.centroid = UBS93.data$UBS93.centroid,
                                                          marker.gene = UBS93.data$UBS93.gene.df$ENSEMBL,
                                                          HER2.amp.signature.genes = UBS93.data$HER2.amp.signature.genes,
                                                          ssgsea.cutoff = UBS93.data$ssgsea.cutoff)

HCC70_SHELF3_pseudobulk_cor               <- predictor.res$cor.matrix %>% as.data.frame()
HCC70_SHELF3_pseudobulk_cor_cpm           <- cbind(HCC70_SHELF3_pseudobulk_cor,HCC70_SHELF3_count_sum_log2cpm["ENSG00000163435",])
colnames(HCC70_SHELF3_pseudobulk_cor_cpm)[4]  <- "ELF3_pseudobulk_log2cpm"
saveRDS(HCC70_SHELF3_pseudobulk_cor_cpm,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_SHELF3_pseudobulk_cor_cpm.rds")
library(ggpubr)
ggplot(HCC70_SHELF3_pseudobulk_cor_cpm,aes(x=ELF3_pseudobulk_log2cpm,y=Basal))+geom_point()+
  geom_smooth(method = "lm", color = "red",se = F)+
  stat_cor(method = "spearman", label.x = 5, label.y = 0.4,size=5)
ggplot(HCC70_SHELF3_pseudobulk_cor_cpm,aes(x=ELF3_pseudobulk_log2cpm,y=Claudin_low))+geom_point()+
  geom_smooth(method = "lm", color = "red",se = F)+
  stat_cor(method = "spearman", label.x = 5, label.y = -0.1,size=5)


HCC70_SHELF3_count_sum_log2cpm_df  <- t(HCC70_SHELF3_count_sum_log2cpm)  %>% as.data.frame()
p1 <- ggplot(HCC70_SHELF3_count_sum_log2cpm_df,aes(x=ENSG00000163435,y=ENSG00000165215))+geom_point(size=1)+
  geom_smooth(method = "lm", color = "red",se = F)+
  stat_cor(method = "spearman", label.x = 5, label.y = 0.4,size=5)+xlab("ELF3")+ylab("CLDN3")
p2 <- ggplot(HCC70_SHELF3_count_sum_log2cpm_df,aes(x=ENSG00000163435,y=ENSG00000189143))+geom_point(size=1)+
  geom_smooth(method = "lm", color = "red",se = F)+
  stat_cor(method = "spearman", label.x = 5, label.y = -0.1,size=5)+xlab("ELF3")+ylab("CLDN4")
p3 <- ggplot(HCC70_SHELF3_count_sum_log2cpm_df,aes(x=ENSG00000163435,y=ENSG00000181885))+geom_point(size=1)+
  geom_smooth(method = "lm", color = "red",se = F)+
  stat_cor(method = "spearman", label.x = 5, label.y = -0.1,size=5)+xlab("ELF3")+ylab("CLDN7")
p4 <- ggplot(HCC70_SHELF3_count_sum_log2cpm_df,aes(x=ENSG00000163435,y=ENSG00000119888))+geom_point(size=1)+
  geom_smooth(method = "lm", color = "red",se = F)+
  stat_cor(method = "spearman", label.x = 5, label.y = -0.1,size=5)+xlab("ELF3")+ylab("EPCAM")
p5 <- ggplot(HCC70_SHELF3_count_sum_log2cpm_df,aes(x=ENSG00000163435,y=ENSG00000026025))+geom_point(size=1)+
  geom_smooth(method = "lm", color = "red",se = F)+
  stat_cor(method = "spearman", label.x = 5, label.y = -0.1,size=5)+xlab("ELF3")+ylab("VIM")
p6 <- ggplot(HCC70_SHELF3_count_sum_log2cpm_df,aes(x=ENSG00000163435,y=ENSG00000186847))+geom_point(size=1)+
  geom_smooth(method = "lm", color = "red",se = F)+
  stat_cor(method = "spearman", label.x = 5, label.y = -0.1,size=5)+xlab("ELF3")+ylab("KRT14")
p7 <- ggplot(HCC70_SHELF3_count_sum_log2cpm_df,aes(x=ENSG00000163435,y=ENSG00000186081))+geom_point(size=1)+
  geom_smooth(method = "lm", color = "red",se = F)+
  stat_cor(method = "spearman", label.x = 5, label.y = -0.1,size=5)+xlab("ELF3")+ylab("KRT5")

combined_plot       <- plot_grid(p1, p2,p3,p4,p5,p6,p7,ncol = 3) 
combined_plot


ggplot(HCC70_SHELF3_cpm,aes(x=ENSG00000163435,y=ENSG00000186847))+geom_point()+
  geom_smooth(method = "lm", color = "red",se = F)+
  stat_cor(method = "spearman", label.x = 0.1, label.y = 3,size=5)+xlab("ELF3")+ylab("KRT14")
ggplot(HCC70_SHELF3_cpm,aes(x=ENSG00000163435,y=ENSG00000186081))+geom_point()+
  geom_smooth(method = "lm", color = "red",se = F)+
  stat_cor(method = "spearman", label.x = 0.1, label.y = 3.5,size=5)+xlab("ELF3")+ylab("KRT5")

#### calculate the correlation coefficients between Claudinn_low and CCLE cell lines####
Claudin_low_scRNA             <- subset(HCC70_SHELF3_scRNA,subtype=="Claudin_low") 
Idents(Claudin_low_scRNA)     <- "orig.ident"
Claudin_low_pesudobulk        <- AggregateExpression(
  Claudin_low_scRNA,
  assays = NULL,
  features = NULL,
  return.seurat = FALSE,
  group.by = "ident",
  add.ident = NULL,
  slot = "count",
  verbose = TRUE
)
Claudin_low_pesudobulk         <- Claudin_low_pesudobulk$RNA
library(edgeR)
Claudin_low_pesudobulk_cpm     <- cpm(Claudin_low_pesudobulk,log=T)

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
CCLE.rna.seq.marker.gene.1000                 <- names(sort(rank.sd,decreasing =TRUE))[1:1000]



####* Function to pick out cell line ##########
library(foreach)
pick.out.cell.line <- function(expr.of.samples,expr.of.cell.lines,marker.gene){
  marker.gene           <- intersect(rownames(expr.of.samples),(marker.gene))
  marker.gene           <- intersect(rownames(expr.of.cell.lines),(marker.gene))
  correlation.matrix    <- cor(expr.of.samples[marker.gene,],expr.of.cell.lines[marker.gene,],method='spearman')
  cell.line.median.cor  <- apply(correlation.matrix,2,function(x) median(x[is.na(x) == FALSE])) %>% sort(decreasing = TRUE)
  best.cell.line        <- names(cell.line.median.cor)[1]
  p.value.vec           <- foreach(cell.line= setdiff(names(cell.line.median.cor),best.cell.line),.combine='c') %do% { #Find the different elements in vector x and vector y (only the different elements in x are taken)
    v                     <- correlation.matrix[,cell.line]
    p.value               <- wilcox.test(correlation.matrix[,best.cell.line],v,alternative = 'greater',paired = TRUE)$p.value
  }
  names(p.value.vec)    <- setdiff(names(cell.line.median.cor),best.cell.line)
  fdr.vec               <- p.adjust(p.value.vec,method='fdr')
  list(cell.line.median.cor=cell.line.median.cor,best.cell.line=best.cell.line,compare.fdr.vec=fdr.vec,correlation.matrix = correlation.matrix )
}


rs.Claudin.low                               <- pick.out.cell.line(expr.of.samples = Claudin_low_pesudobulk_cpm,expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor                                     <- rs.Claudin.low$correlation.matrix  %>% t() %>% as.data.frame()
colnames(expr.cor)[1]  <- "cor"
expr.cor$CCLE_Name     <- rownames(expr.cor)
brca_info              <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_expr.cor          <- expr.cor[expr.cor$CCLE_name %in% brca_meta$cell.line.name,]
VlnPlot(Claudin_low_scRNA,features = c("nCount_RNA"),pt.size = 3)
write.csv(brca_expr.cor,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/brca_expr.cor.csv")
write.csv(expr.cor,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/expr.cor.csv")


rs.Claudin.low.scdata                               <- pick.out.cell.line(expr.of.samples = as.matrix(Claudin_low_scRNA@assays$RNA@data),expr.of.cell.lines = CCLE.log2.rpkm.matrix,marker.gene = CCLE.rna.seq.marker.gene.1000)
expr.cor.scdata                                     <- rs.Claudin.low.scdata$correlation.matrix  %>% t() %>% as.data.frame()
expr.cor.scdata$median_cor    <- apply(expr.cor.scdata,1,median)
expr.cor.scdata$CCLE_Name     <- rownames(expr.cor.scdata)
brca_info                     <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca.expr.cor.scdata          <- merge(expr.cor.scdata,brca_info[,c(2,22)],by="CCLE_Name")



#########
df                            <- data.frame(logELF3=as.vector(HCC70_SHELF3_scRNA@assays$RNA@data["ENSG00000163435",]),log.count=log1p(HCC70_SHELF3_scRNA$nCount_RNA))
ggplot(df,aes(x=log.count,y=logELF3))+geom_point()+
  geom_smooth(method = "loess", color = "red",se = F)+ggtitle("HCC70_SHELF3")

df_SCR                            <- data.frame(logELF3=as.vector(HCC70_SCR_scRNA@assays$RNA@data["ENSG00000163435",]),log.count=log1p(HCC70_SCR_scRNA$nCount_RNA))
ggplot(df_SCR,aes(x=log.count,y=logELF3))+geom_point()+ggtitle("control_HCC70")



####DE analysis for after filtering the cells having the low count####
HCC70_SHELF3_scRNA$log1pCount    <- log1p(HCC70_SHELF3_scRNA$nCount_RNA)
HCC70_ELF3_highCount_scRNA       <- subset(HCC70_SHELF3_scRNA,log1pCount >10 )
HCC70_SHELF3_cpm                 <- HCC70_ELF3_highCount_scRNA@assays$RNA@data %>% as.matrix()%>% t() %>% as.data.frame()
HCC70_SHELF3_cpm                 <- HCC70_SHELF3_cpm[order(HCC70_SHELF3_cpm$ENSG00000163435,decreasing = F),]
all(diff(HCC70_SHELF3_cpm[,"ENSG00000163435"]) >= 0)
HCC70_SHELF3_cpm_high_low        <- HCC70_SHELF3_cpm[c(1:100,6573:6672),]
HCC70_SHELF3_high_low_scRNA      <- HCC70_SHELF3_scRNA[,rownames(HCC70_SHELF3_cpm_high_low)]
HCC70_SHELF3_high_low_scRNA$group    <- c(rep("ELF3_low",100),rep("ELF3_high",100))


Idents(HCC70_SHELF3_high_low_scRNA)  <- "group"
#wilcox
HCC70_SHELF3_high_low_de_res            <- FindMarkers(HCC70_SHELF3_high_low_scRNA,slot = "data",test.use = "wilcox",ident.1="ELF3_low",ident.2="ELF3_high",logfc.threshold = 0,pseudocount.use = 10^(-9))
HCC70_SHELF3_high_low_deres             <- within(HCC70_SHELF3_high_low_de_res,{
  sig                                   <- NA
  sig[p_val_adj < 0.05 & avg_log2FC > 0]   <- "up"
  sig[p_val_adj < 0.05 & avg_log2FC < 0]  <- "down"
})
HCC70_SHELF3_high_low_deres[is.na(HCC70_SHELF3_high_low_deres$sig),6] <- "none"
#
symbol                                <- bitr(rownames(HCC70_SHELF3_high_low_deres),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
HCC70_SHELF3_high_low_deres$ENSEMBL   <- rownames(HCC70_SHELF3_high_low_deres)
HCC70_SHELF3_high_low_deres_symbol    <- merge(HCC70_SHELF3_high_low_deres,symbol,by="ENSEMBL")
HCC70_SHELF3_low_VS_high_deres_symbol_highCount <- HCC70_SHELF3_high_low_deres_symbol
write.csv(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_SHELF3_low_VS_high_deres_symbol_highCount.csv")
saveRDS(HCC70_SHELF3_high_low_scRNA,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_SHELF3_high_low_scRNA.rds")

HCC70_SHELF3_high_low_exp_df        <- HCC70_SHELF3_high_low_scRNA@assays$RNA@data %>% as.matrix() %>% t()%>% as.data.frame()
HCC70_SHELF3_high_low_exp_df$group  <- c(rep("ELF3_low",100),rep("ELF3_high",100))
p1  <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000165215))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("CLDN3")
p2 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000189143))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("CLDN4")
p3 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000181885))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("CLDN7")
p4 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000119888))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("EPCAM")
p5 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000026025))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("VIM")
p6 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000163435))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("ELF3")
p7 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000186847))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("KRT14")
p8 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=ENSG00000186081))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("KRT5")
combined_plot       <- plot_grid(p1, p2,p3 ,p4,p5,p6,p7,p8,ncol = 3) 
combined_plot


HCC70_SHELF3_high_low_exp_df$Basal  <- HCC70_SHELF3_high_low_scRNA$Basal
HCC70_SHELF3_high_low_exp_df$Claudin_low <- HCC70_SHELF3_high_low_scRNA$Claudin_low
HCC70_SHELF3_high_low_exp_df$HORL   <- HCC70_SHELF3_high_low_scRNA$H_or_L
p1 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=Basal))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("Basal")
p2 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=Claudin_low))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("Claudin_low")
p3 <- ggplot(HCC70_SHELF3_high_low_exp_df,aes(x=group,y=HORL))+
  geom_boxplot() +stat_compare_means(method = "wilcox.test")+ylab("HORL")
combined_plot       <- plot_grid(p1, p2,p3 ,ncol = 3) 
combined_plot




####intergrate HCC70 SCR and HCC70 SHELF3####
HCC70_SCR_scRNA$orig.ident     <- rep("HCC70_SCR_scRNA",ncol(HCC70_SCR_scRNA))
rown                           <- paste0("HCC70_SCR_",colnames(HCC70_SCR_scRNA))
names(rown)                    <- Cells(HCC70_SCR_scRNA)
HCC70_SCR_scRNA                <- RenameCells(HCC70_SCR_scRNA, new.names = rown)

HCC70_SHELF3_scRNA$orig.ident    <- rep("HCC70_SHELF3_scRNA",ncol(HCC70_SHELF3_scRNA))
rown_1                           <- paste0("HCC70_SHELF3_",colnames(HCC70_SHELF3_scRNA))

names(rown_1)                    <- Cells(HCC70_SHELF3_scRNA)
HCC70_SHELF3_scRNA                <- RenameCells(HCC70_SHELF3_scRNA, new.names = rown_1)

HCC70_scRNA                    <- merge(HCC70_SHELF3_scRNA,HCC70_SCR_scRNA)

#features choosing
HCC70_scRNA                     <- FindVariableFeatures(HCC70_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                               <- head(VariableFeatures(HCC70_scRNA),10)

#data scaling 
all_gene                               <- rownames(HCC70_scRNA)
HCC70_scRNA                     <- ScaleData(HCC70_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data

#linear dimensionality reduction
HCC70_scRNA                     <- RunPCA(HCC70_scRNA,features = VariableFeatures(object=HCC70_scRNA))

#dimension selection
HCC70_scRNA                     <- JackStraw(HCC70_scRNA,num.replicate = 100)
HCC70_scRNA                     <- ScoreJackStraw(HCC70_scRNA,dims = 1:20)
JackStrawPlot(HCC70_scRNA,dims=1:20)
ElbowPlot(HCC70_scRNA,ndims = 50)#choose dims of 30
#cell clustering
HCC70_scRNA                     <- FindNeighbors(HCC70_scRNA,dims = 1:40)
HCC70_scRNA                     <- FindClusters(HCC70_scRNA,resolution = 0.1)
table(Idents(HCC70_scRNA))

#unlinear dimension reduction(UMAP or tSNE)
HCC70_scRNA                     <- RunUMAP(HCC70_scRNA,dims = 1:40)
HCC70_scRNA                     <- RunTSNE(HCC70_scRNA,dims = 1:40,check_duplicates = FALSE)
Idents(HCC70_scRNA)             <- "orig.ident"
DimPlot(HCC70_scRNA,reduction = "umap",label = T)
DimPlot(HCC70_scRNA,reduction = "tsne",label = T)
saveRDS(HCC70_scRNA,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_scRNA.rds")
DimPlot(HCC70_scRNA,reduction = "pca",label = T)



####Explore whether it is due to low count that those parts are clustered together.####
HCC70_SHELF3_cluster1_scRNA                      <- subset(HCC70_SHELF3_scRNA,seurat_clusters=="1")
set.seed(456)
HCC70_SHELF3_cluster1_200cells                   <- HCC70_SHELF3_cluster1_scRNA[,sample(1:ncol(HCC70_SHELF3_cluster1_scRNA),200,replace = F)]
boxplot(HCC70_SHELF3_cluster1_200cells$log1pCount)
median(HCC70_SHELF3_cluster1_200cells$log1pCount)

library(scuttle)
HCC70_SHELF3_cluster1_200cells_count             <- HCC70_SHELF3_cluster1_200cells@assays$RNA@counts  %>% as.matrix() 
HCC70_SHELF3_cluster1_200cells_count_downsampled <- downsampleMatrix(HCC70_SHELF3_cluster1_200cells_count,prop = 0.06)


HCC70_downsampled_scRNA                    <- CreateSeuratObject(counts = HCC70_SHELF3_cluster1_200cells_count_downsampled,min.cells = 3, min.features = 200)
#calculate the proportion of ribosomal genes in cells
HCC70_downsampled_scRNA[["percent.mt"]]     <- PercentageFeatureSet(HCC70_downsampled_scRNA,pattern = "^MT-")
head(HCC70_downsampled_scRNA@meta.data,5)
VlnPlot(HCC70_downsampled_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
HCC70_downsampled_scRNA                     <- subset(HCC70_downsampled_scRNA,nFeature_RNA >0 )

#normalized
HCC70_downsampled_scRNA                     <- NormalizeData(HCC70_downsampled_scRNA)

#merge data
HCC70_downsampled_scRNA$log1pCount   <- log1p(HCC70_downsampled_scRNA$nCount_RNA)
HCC70_downsampled_scRNA$orig.ident   <- rep("downsampled",ncol(HCC70_downsampled_scRNA))
rown_1                               <- paste0("Down",colnames(HCC70_downsampled_scRNA))
names(rown_1)                        <- Cells(HCC70_downsampled_scRNA)
HCC70_downsampled_scRNA              <- RenameCells(HCC70_downsampled_scRNA, new.names = rown_1)
HCC70_SHELF3_integrated_scRNA        <- merge(HCC70_SHELF3_scRNA,HCC70_downsampled_scRNA)
#features choosing
HCC70_SHELF3_integrated_scRNA              <- FindVariableFeatures(HCC70_SHELF3_integrated_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.

#data scaling 
all_gene                                   <- rownames(HCC70_SHELF3_integrated_scRNA)
HCC70_SHELF3_integrated_scRNA              <- ScaleData(HCC70_SHELF3_integrated_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data

#linear dimensionality reduction
HCC70_SHELF3_integrated_scRNA              <- RunPCA(HCC70_SHELF3_integrated_scRNA,features = VariableFeatures(object=HCC70_SHELF3_integrated_scRNA))

#dimension selection
HCC70_SHELF3_integrated_scRNA              <- JackStraw(HCC70_SHELF3_integrated_scRNA,num.replicate = 100)
HCC70_SHELF3_integrated_scRNA              <- ScoreJackStraw(HCC70_SHELF3_integrated_scRNA,dims = 1:20)
JackStrawPlot(HCC70_SHELF3_integrated_scRNA,dims=1:20)
ElbowPlot(HCC70_SHELF3_integrated_scRNA,ndims = 50)#choose dims of 30
#cell clustering
HCC70_SHELF3_integrated_scRNA                     <- FindNeighbors(HCC70_SHELF3_integrated_scRNA,dims = 1:40)
HCC70_SHELF3_integrated_scRNA                     <- FindClusters(HCC70_SHELF3_integrated_scRNA,resolution = 0.1)
table(Idents(HCC70_SHELF3_integrated_scRNA))

#unlinear dimension reduction(UMAP or tSNE)
HCC70_SHELF3_integrated_scRNA                     <- RunUMAP(HCC70_SHELF3_integrated_scRNA,dims = 1:40)
Idents(HCC70_SHELF3_integrated_scRNA)             <- "orig.ident"
DimPlot(HCC70_SHELF3_integrated_scRNA,reduction = "umap",label = T)

####*TSNE####
HCC70_SHELF3_integrated_exp                       <- as.matrix(HCC70_SHELF3_integrated_scRNA@assays$RNA@data)
HCC70_SHELF3_integrated_exp_UBS83                 <- HCC70_SHELF3_integrated_exp[rownames(HCC70_SHELF3_integrated_exp) %in% UBS93.data$UBS93.gene.df$ENSEMBL,]
library(Rtsne)
cor_matrix                                        <- cor(HCC70_SHELF3_integrated_exp_UBS83,method = "spearman")
dissimilarity_matrix                              <- 1 - cor_matrix
normalized_dissimilarity_matrix                   <- dissimilarity_matrix / max(dissimilarity_matrix)
HCC70_SHELF3_integrated_scRNA                     <- RunTSNE(
  HCC70_SHELF3_integrated_scRNA,
  reduction = "pca",
  cells = NULL,
  dims = 1:40,
  features = NULL,
  seed.use = 1,
  tsne.method = "Rtsne",
  dim.embed = 2,
  distance.matrix = normalized_dissimilarity_matrix,
  reduction.name = "tsne",
  reduction.key = "tSNE_"
)
DimPlot(HCC70_SHELF3_integrated_scRNA,reduction = "tsne",label = T)


cluster_2 <- subset(HCC70_SHELF3_scRNA,seurat_clusters=="2")
cluster_2 <- colnames(cluster_2)
df        <- data.frame(id=colnames(HCC70_SHELF3_scRNA),ident=HCC70_SHELF3_scRNA$seurat_clusters)
df_2      <- data.frame(id=colnames(HCC70_downsampled_scRNA),ident=HCC70_downsampled_scRNA$orig.ident)
df        <- rbind(df,df_2)
identical(df$id,colnames(HCC70_SHELF3_integrated_scRNA))
HCC70_SHELF3_integrated_scRNA$cluster  <- df$ident
Idents(HCC70_SHELF3_integrated_scRNA)  <- "cluster"
DimPlot(HCC70_SHELF3_integrated_scRNA,reduction = "umap",label = T,pt.size = 2)
DimPlot(HCC70_SHELF3_integrated_scRNA,reduction = "tsne",label = T,pt.size = 3)
saveRDS(HCC70_SHELF3_integrated_scRNA,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_SHELF3_integrated_scRNA.rds")


library(Seurat)
expr_matrix     <- HCC70_SHELF3_scRNA[["RNA"]]@counts
zero_counts     <- apply(expr_matrix == 0, 2, sum)
zero_ratio      <- zero_counts / nrow(expr_matrix)
HCC70_SHELF3_scRNA$zero_ratio  <- zero_ratio

expr_matrix     <- HCC70_SHELF3_integrated_scRNA[["RNA"]]@counts
zero_counts     <- apply(expr_matrix == 0, 2, sum)
zero_ratio      <- zero_counts / nrow(expr_matrix)
HCC70_SHELF3_integrated_scRNA$zero_ratio  <- zero_ratio

saveRDS(HCC70_SHELF3_scRNA,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_SHELF3_scRNA.rds")


####down genes in CCLE ####
load("Pro_TNBC/data/CCLE/CCLE.RData")
brca_info                  <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_info                  <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
down_gene_id               <- subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,sig=="down")
brca_log2rpkm              <- CCLE.log2.rpkm.matrix[,colnames(CCLE.log2.rpkm.matrix) %in% brca_info$CCLE_Name]
brca_log2rpkm_down         <- brca_log2rpkm[rownames(brca_log2rpkm) %in% down_gene_id$ENSEMBL,] %>% as.data.frame()
brca_log2rpkm_down$ENSEMBL <- rownames(brca_log2rpkm_down)
brca_log2rpkm_down         <- merge(down_gene_id[,c(2,9)],brca_log2rpkm_down,by="ENSEMBL")
brca_log2rpkm_down         <- column_to_rownames(brca_log2rpkm_down,var = "SYMBOL")
brca_log2rpkm_down         <- brca_log2rpkm_down[,-1] %>% as.matrix() %>% t()
brca_info                  <- arrange(brca_info,lineage_molecular_subtype)
brca_log2rpkm_down                   <- brca_log2rpkm_down[match(brca_info$CCLE_Name,rownames(brca_log2rpkm_down)),]
identical(rownames(brca_log2rpkm_down),brca_info$CCLE_Name)
annotation.col              <- brca_info[,22,drop=F] %>% as.data.frame()
rownames(annotation.col)    <- brca_info$CCLE_Name 
colnames(annotation.col)[1] <- "subtype"
annot.colors                <- list(subtype=c("basal_A"="Red","basal_B"= "Blue","HER2_amp"= "Purple","luminal"= "#FFCC00","luminal_HER2_amp"="grey"))
a <- names(table(annotation.col))
library(pheatmap)
p1 <- pheatmap(t(as.matrix(brca_log2rpkm_down)),cluster_rows = T,
               cluster_cols = F,
               color = colorRampPalette(c("navy","white","firebrick3"))(100),
               show_colnames = F,border_color = NA,
               scale = "row",
               show_rownames = T,annotation_col   = annotation.col,
               annotation_colors = annot.colors,fontsize=30,cellwidth = 18, cellheight = 25,
               gaps_col =c(13,28,37,48)
)
ggsave(p1,filename = "Pro_TNBC/paper/plot/section_5/heatmap.pdf",width = 25,height = 40)




####up genes in CCLE ####
load("Pro_TNBC/data/CCLE/CCLE.RData")
brca_info                  <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_info                  <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
up_gene_id               <- subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,sig=="up")
brca_log2rpkm              <- CCLE.log2.rpkm.matrix[,colnames(CCLE.log2.rpkm.matrix) %in% brca_info$CCLE_Name]
brca_log2rpkm_up         <- brca_log2rpkm[rownames(brca_log2rpkm) %in% up_gene_id$ENSEMBL,] %>% as.data.frame()
brca_log2rpkm_up$ENSEMBL <- rownames(brca_log2rpkm_up)
brca_log2rpkm_up         <- merge(up_gene_id[,c(2,9)],brca_log2rpkm_up,by="ENSEMBL")
brca_log2rpkm_up         <- column_to_rownames(brca_log2rpkm_up,var = "SYMBOL")
brca_log2rpkm_up         <- brca_log2rpkm_up[,-1] %>% as.matrix() %>% t()
brca_info                  <- arrange(brca_info,lineage_molecular_subtype)
brca_log2rpkm_up                   <- brca_log2rpkm_up[match(brca_info$CCLE_Name,rownames(brca_log2rpkm_up)),]
identical(rownames(brca_log2rpkm_up),brca_info$CCLE_Name)
annotation.col              <- brca_info[,22,drop=F] %>% as.data.frame()
rownames(annotation.col)    <- brca_info$CCLE_Name 
colnames(annotation.col)[1] <- "subtype"
annot.colors                <- list(subtype=c("basal_A"="Red","basal_B"= "Blue","HER2_amp"= "Purple","luminal"= "#FFCC00","luminal_HER2_amp"="grey"))
library(pheatmap)
p1 <- pheatmap(t(as.matrix(brca_log2rpkm_up)),cluster_rows = T,
               cluster_cols = F,
               color = colorRampPalette(c("navy","white","firebrick3"))(100),
               show_colnames = F,border_color = NA,
               scale = "row",
               show_rownames = T,annotation_col   = annotation.col,
               annotation_colors = annot.colors,fontsize=30,cellwidth = 18, cellheight = 25,
               gaps_col =c(13,28,37,48)
)
ggsave(p1,filename = "Pro_TNBC/paper/plot/section_5/heatmap_up.pdf",width = 25,height = 100,limitsize = F)


brca_basal_cl               <- subset(brca_info,lineage_molecular_subtype=="basal_A"|lineage_molecular_subtype=="basal_B")
brca_log2rpkm_up_bc         <- brca_log2rpkm_up[rownames(brca_log2rpkm_up) %in% brca_basal_cl$CCLE_Name,]
identical(rownames(brca_log2rpkm_up_bc),brca_basal_cl$CCLE_Name)
annotation.col              <- brca_basal_cl[,22,drop=F] %>% as.data.frame()
rownames(annotation.col)    <- brca_basal_cl$CCLE_Name 
colnames(annotation.col)[1] <- "subtype"
annot.colors                <- list(subtype=c("basal_A"="Red","basal_B"= "Blue"))
library(pheatmap)
p2 <- pheatmap(t(as.matrix(brca_log2rpkm_up_bc)),cluster_rows = T,
               cluster_cols = F,
               color = colorRampPalette(c("navy","white","firebrick3"))(100),
               show_colnames = F,border_color = NA,
               scale = "row",
               show_rownames = T,annotation_col   = annotation.col,
               annotation_colors = annot.colors,fontsize=30,cellwidth = 18, cellheight = 25
)
ggsave(p2,filename = "Pro_TNBC/paper/plot/section_5/heatmap_up_bc.pdf",width = 25,height = 80,limitsize = F)



tab
p1 <- ggplot(brca_log2rpkm_down_bc_df,aes(x=subtype,y=KRT19))+geom_boxplot()+geom_point(position = "jitter")+stat_compare_means(method = "wilcox.test")
p2 <- ggplot(brca_log2rpkm_down_bc_df,aes(x=subtype,y=CLDN4))+geom_boxplot()+geom_point(position = "jitter")+stat_compare_means(method = "wilcox.test")
p3 <- ggplot(brca_log2rpkm_down_bc_df,aes(x=subtype,y=ELF3))+geom_boxplot()+geom_point(position = "jitter")+stat_compare_means(method = "wilcox.test")
p4 <- ggplot(brca_log2rpkm_down_bc_df,aes(x=subtype,y=RHOV))+geom_boxplot()+geom_point(position = "jitter")+stat_compare_means(method = "wilcox.test")
p5 <- ggplot(brca_log2rpkm_down_bc_df,aes(x=subtype,y=TMPRSS2))+geom_boxplot()+geom_point(position = "jitter")+stat_compare_means(method = "wilcox.test")
p6 <- ggplot(brca_log2rpkm_down_bc_df,aes(x=subtype,y=LOC124905021))+geom_boxplot()+geom_point(position = "jitter")+stat_compare_means(method = "wilcox.test")
p7 <- ggplot(brca_log2rpkm_down_bc_df,aes(x=subtype,y=SPNS2))+geom_boxplot()+geom_point(position = "jitter")+stat_compare_means(method = "wilcox.test")
p8 <- ggplot(brca_log2rpkm_down_bc_df,aes(x=subtype,y=TACSTD2))+geom_boxplot()+geom_point(position = "jitter")+stat_compare_means(method = "wilcox.test")
p9 <- ggplot(brca_log2rpkm_down_bc_df,aes(x=subtype,y=VTCN1))+geom_boxplot()+geom_point(position = "jitter")+stat_compare_means(method = "wilcox.test")

combined_plot       <- plot_grid(p1, p2,p3,p4,p5,p6,p7,p8,p9, ncol = 3) 
combined_plot


brca_log2rpkm_exp_bc            <- brca_log2rpkm[,colnames(brca_log2rpkm) %in% brca_basal_cl$CCLE_Name] %>% t()
brca_log2rpkm_exp_bc_df         <- brca_log2rpkm_exp_bc[match(brca_basal_cl$CCLE_Name,rownames(brca_log2rpkm_exp_bc)),] %>% as.data.frame()
identical(rownames(brca_log2rpkm_down_bc_df),brca_basal_cl$CCLE_Name)
brca_log2rpkm_exp_bc_df$subtype <- brca_basal_cl$lineage_molecular_subtype
ggplot(brca_log2rpkm_exp_bc_df,aes(x=subtype,y=ENSG00000026025))+geom_boxplot()+geom_point(position = "jitter")+stat_compare_means(method = "wilcox.test")

#### AUC score####
# AUCell score ------------------------------------------------------------


basal_marker              <- c("KRT5","KRT14","TP63","ACTA2","MYLK")
ENSEMBL                   <- bitr(basal_marker,fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = "org.Hs.eg.db")
basal_marker              <- ENSEMBL$ENSEMBL
sig                       <- list(basal_marker=basal_marker)
cell_rangkings                                                          <- AUCell_buildRankings(HCC70_SHELF3_high_low_scRNA[['RNA']]@data,nCores = 30)
cells_AUC                                                               <- AUCell_calcAUC(msigdbr.geneset.list, cell_rangkings,nCores = 5)


auc_score                                  <- t(getAUC(cells_AUC))
auc_score_df                               <- as.data.frame(auc_score)
auc_score_df$group                         <- HCC70_SHELF3_high_low_scRNA$group

p.test <- foreach(i = colnames(auc_score_df[,1:50]),.combine=c)%do%{
  value   <- colnames(auc_score_df[,1:50])[i]
  res  <- wilcox.test(value~ group,data = auc_score_df)
  p    <- res$p.value
}
HCC70_SHELF3_high_low_scRNA                      <- AddMetaData(HCC70_SHELF3_high_low_scRNA, auc_score_df)
ggplot(HCC70_SHELF3_high_low_scRNA@meta.data,aes(x=group,y=))+geom_boxplot()
BiocManager::install("msigdbr")
library(msigdbr)
msigdbr.geneset  <- msigdbr(category = 'H')
msigdbr.geneset.list  <-split(msigdbr.geneset$ensembl_gene,msigdbr.geneset$gs_name) 

wilcox_results  <- lapply(1:50, function(i) {
  value_col     <- auc_score_df[, i]
  test_result   <- wilcox.test(value_col ~ group, data = auc_score_df)
  p_value = test_result$p.value
})
wilcox_results               <- unlist(wilcox_results) %>% as.data.frame()
wilcox_results$pathway       <- colnames(auc_score_df)[1:50]
colnames(wilcox_results)[1]  <- "p.value"
wilcox_results$fdr           <- p.adjust(wilcox_results$p.value,method = "fdr")
wilcox_results_de            <- subset(wilcox_results,fdr < 0.05)
auc_score_df_de <- auc_score_df[,colnames(auc_score_df)%in% c(wilcox_results_de$pathway,"group")]
pdf("Pro_TNBC/paper/plot/section_5/pathway.boxplots.pdf")
lapply(1:38, function(i) {
  value_col_name <- colnames(auc_score_df_de)[i]
  ggplot(auc_score_df_de,aes(y=auc_score_df_de[, i] ,x=group))+
    geom_boxplot()+stat_compare_means(method = "wilcox.test")+ggtitle(value_col_name) 
  
})
dev.off()





TNBC_cell_rangkings                                                          <- AUCell_buildRankings(TNBC_BRCA1_0554_tumor_scRNA_BC[['RNA']]@data,nCores = 30)
TNBC_cells_AUC                                                               <- AUCell_calcAUC(msigdbr.geneset.list, TNBC_cell_rangkings,nCores = 5)


TNBC_auc_score                                  <- t(getAUC(TNBC_cells_AUC))
TNBC_auc_score_df                               <- as.data.frame(TNBC_auc_score)
TNBC_auc_score_df$subtype                         <- TNBC_BRCA1_0554_tumor_scRNA_BC$subtype


TNBC_wilcox_results  <- lapply(1:50, function(i) {
  value_col     <- TNBC_auc_score_df[, i]
  test_result   <- wilcox.test(value_col ~ subtype, data = TNBC_auc_score_df)
  p_value = test_result$p.value
})
TNBC_wilcox_results               <- unlist(TNBC_wilcox_results) %>% as.data.frame()
TNBC_wilcox_results$pathway       <- colnames(TNBC_auc_score_df)[1:50]
colnames(TNBC_wilcox_results)[1]  <- "p.value"
TNBC_wilcox_results$fdr           <- p.adjust(TNBC_wilcox_results$p.value,method = "fdr")
TNBC_wilcox_results_de            <- subset(TNBC_wilcox_results,fdr < 0.05)
TNBC_auc_score_df_de <- TNBC_auc_score_df[,colnames(TNBC_auc_score_df)%in% c(TNBC_wilcox_results_de$pathway,"subtype")]
pdf("Pro_TNBC/paper/plot/section_5/TNBC0554.pathway.boxplots.pdf")
lapply(1:20, function(i) {
  value_col_name <- colnames(TNBC_auc_score_df_de)[i]
  ggplot(TNBC_auc_score_df_de,aes(y=TNBC_auc_score_df_de[, i] ,x=subtype))+
    geom_boxplot()+stat_compare_means(method = "wilcox.test")+
    ggtitle(value_col_name) 
  
})
dev.off()

intersect_pathway <- intersect(TNBC_wilcox_results_de$pathway,wilcox_results_de$pathway)



####DE analysis for SCR and SHELF3 group ####
Idents(HCC70_scRNA)  <- "orig.ident"
#wilcox
HCC70_de_res            <- FindMarkers(HCC70_scRNA,slot = "data",test.use = "wilcox",ident.1="HCC70_SHELF3_scRNA",ident.2="HCC70_SCR_scRNA",logfc.threshold = 0,pseudocount.use = 10^(-9))
HCC70_deres             <- within(HCC70_de_res,{
  sig                                   <- NA
  sig[p_val_adj < 0.05 & avg_log2FC > 1]   <- "up"
  sig[p_val_adj < 0.05 & avg_log2FC < -1]  <- "down"
})
HCC70_deres[is.na(HCC70_deres$sig),6] <- "none"

symbol                                <- bitr(rownames(HCC70_deres),fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
HCC70_deres$ENSEMBL   <- rownames(HCC70_deres)
HCC70_deres_symbol    <- merge(HCC70_deres,symbol,by="ENSEMBL")
HCC70_SHELF3_VS_SCR_deres_symbol <- HCC70_deres_symbol
write.csv(HCC70_SHELF3_VS_SCR_deres_symbol,file = "Pro_TNBC/paper/data/results/section_5/ourPipeline_HCC70/HCC70_SHELF3_VS_SCR_deres_symbol.csv")





HCC70_SHELF3_low_cpm                   <- subset(HCC70_SHELF3_high_low_scRNA,group=="ELF3_low")
HCC70_SHELF3_low_cpm                   <- as.matrix(HCC70_SHELF3_low_cpm@assays$RNA@data)
HCC70_SHELF3_low_cpm_scale             <- apply(HCC70_SHELF3_low_cpm, 1,scale)
HCC70_SHELF3_low_cpm_scale_median      <- apply(HCC70_SHELF3_low_cpm_scale,2,median)
HCC70_SHELF3_low_cpm_scale_median_nona <- HCC70_SHELF3_low_cpm_scale_median[!is.na(HCC70_SHELF3_low_cpm_scale_median)]
plot(sort(HCC70_SHELF3_low_cpm_scale_median_nona))
hist(HCC70_SHELF3_low_cpm_scale_median_nona,breaks = 40)
HCC70_SHELF3_low_cpm_scale_median_nona <- HCC70_SHELF3_low_cpm_scale_median_nona[sort(HCC70_SHELF3_low_cpm_scale_median_nona)]
all(diff(HCC70_SHELF3_low_cpm_scale_median_nona >=0))
cl_low_genes                           <- names(HCC70_SHELF3_low_cpm_scale_meadian_nona[HCC70_SHELF3_low_cpm_scale_meadian_nona < -0.6])


HCC70_SHELF3_high_cpm                   <- subset(HCC70_SHELF3_high_low_scRNA,group=="ELF3_high")
HCC70_SHELF3_high_cpm                   <- as.matrix(HCC70_SHELF3_high_cpm@assays$RNA@data)
HCC70_SHELF3_high_cpm_scale             <- apply(HCC70_SHELF3_high_cpm, 1,scale)
HCC70_SHELF3_high_cpm_scale_median      <- apply(HCC70_SHELF3_high_cpm_scale,2,median)
HCC70_SHELF3_high_cpm_scale_median_nona <- HCC70_SHELF3_high_cpm_scale_median[!is.na(HCC70_SHELF3_high_cpm_scale_median)]
plot(sort(HCC70_SHELF3_high_cpm_scale_median_nona))
hist(HCC70_SHELF3_high_cpm_scale_median_nona,breaks = 40)

background_genes <- cl_low_genes
group1_genes     <- HCC70_SHELF3_low_VS_high_deres_symbol_highCount[HCC70_SHELF3_low_VS_high_deres_symbol_highCount$sig=="down",]$ENSEMBL 
N                <- length(background_genes)  
K                <- length(intersect(background_genes, group1_genes))  
n                <- length(group1_genes)  
k                <- K  
p_value          <- phyper(k - 1, K, N - K, n, lower.tail = FALSE)
enrichment_score <- -log10(p_value)













####select claudin-low low genes ####
load("Pro_TNBC/data/CCLE/CCLE.RData")
brca_info                             <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_info                             <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
brca_info_cl                          <- subset(brca_info,lineage_molecular_subtype=="basal_B")
brca.log2.rpkm.matrix                 <- CCLE.log2.rpkm.matrix[,colnames(CCLE.log2.rpkm.matrix) %in% brca_info$CCLE_Name]
brca.log2.rpkm.matrix.scale           <- apply(brca.log2.rpkm.matrix,1,scale)
rownames(brca.log2.rpkm.matrix.scale) <- colnames(brca.log2.rpkm.matrix)
brca.log2.rpkm.matrix.CL.scale        <- brca.log2.rpkm.matrix.scale[rownames(brca.log2.rpkm.matrix.scale)%in% brca_info_cl$CCLE_Name,]
brca.log2.rpkm.matrix.CL.scale.median    <- apply(brca.log2.rpkm.matrix.CL.scale,2,median)
brca.log2.rpkm.matrix.CL.scale.median_na <- brca.log2.rpkm.matrix.CL.scale.median[!is.na(brca.log2.rpkm.matrix.CL.scale.median)]
brca.log2.rpkm.matrix.CL.scale.median_na  <- brca.log2.rpkm.matrix.CL.scale.median_na[names(brca.log2.rpkm.matrix.CL.scale.median_na) %in% rownames(HCC70_SHELF3_scRNA)]
plot(sort(brca.log2.rpkm.matrix.CL.scale.median_na))
hist(brca.log2.rpkm.matrix.CL.scale.median_na,breaks = 40)
brca.log2.rpkm.matrix.CL.scale.median_na <- brca.log2.rpkm.matrix.CL.scale.median_na[order(brca.log2.rpkm.matrix.CL.scale.median_na,decreasing = F)]
all(diff(brca.log2.rpkm.matrix.CL.scale.median_na >= 0))
cl_low_genes  <- names(brca.log2.rpkm.matrix.CL.scale.median_na[brca.log2.rpkm.matrix.CL.scale.median_na < -1])

save(cl_low_genes,file = "Pro_TNBC/paper/data/results/section_5/cl_low_genes.RData")


down_genes  <- HCC70_SHELF3_low_VS_high_deres_symbol_highCount[HCC70_SHELF3_low_VS_high_deres_symbol_highCount$sig=="down",]$ENSEMBL 
N           <- length(rownames(HCC70_SHELF3_scRNA))#the number of total genes 
K           <- length(cl_low_genes)#the number of claudin-low cell lines down-regulated genes  
n           <- length(down_genes)#the number of down-regulated in ELF3-low group 
k           <- length(intersect(cl_low_genes, down_genes))  
phyper(q=k-1, m=K, n=N - K, k=n, lower.tail = FALSE)


library(fgsea)
log2fc              <- HCC70_SHELF3_low_VS_high_deres_symbol_highCount
log2fc              <- subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,p_val_adj < 0.05)
rank_vecter         <- log2fc$avg_log2FC
names(rank_vecter)  <- log2fc$ENSEMBL
rank_vecter         <- sort(rank_vecter,decreasing = TRUE)
pathway             <- list(cl_low_genes=cl_low_genes)
fgseaRes            <- fgsea(pathways = pathway, 
                  stats = rank_vecter,
                  eps = 0.0,scoreType="neg")

# plot the most significantly enriched pathway
plotEnrichment(pathway[[1]],rank_vecter) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway)



brca.log2.rpkm.matrix.CL.scale.median_na.df <- as.data.frame(brca.log2.rpkm.matrix.CL.scale.median_na)
brca.log2.rpkm.matrix.CL.scale.median_na.df  <- brca.log2.rpkm.matrix.CL.scale.median_na.df[order(brca.log2.rpkm.matrix.CL.scale.median_na.df$brca.log2.rpkm.matrix.CL.scale.median_na),,drop=F]
brca.log2.rpkm.matrix.CL.scale.median_na.df$rank  <- rank(brca.log2.rpkm.matrix.CL.scale.median_na.df$brca.log2.rpkm.matrix.CL.scale.median_na)
fig1f_1                              <- ggplot(brca.log2.rpkm.matrix.CL.scale.median_na.df,aes(x=rank,y=brca.log2.rpkm.matrix.CL.scale.median_na))+
  geom_point(size=6,show.legend = F)+ggplot.style+xlab(paste("rank",sep = "\n","(n=22713)"))+ylab("Median expression")+
  scale_x_continuous(limits = c(0, 23000),breaks = seq(0,23000,5000))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+geom_hline(aes(yintercept=-1,color="red",linewidth=1))
ggsave(fig1f_1,filename = "Pro_TNBC/paper/plot/method/the.rank.charts.of.clgenes.pdf",width=20,height=15)


####Claudin-low up gene ####
cl_up_genes         <- names(brca.log2.rpkm.matrix.CL.scale.median_na[brca.log2.rpkm.matrix.CL.scale.median_na >1])
cl_up_genes_symbol  <- bitr(cl_up_genes,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")


up_genes    <- HCC70_SHELF3_low_VS_high_deres_symbol_highCount[HCC70_SHELF3_low_VS_high_deres_symbol_highCount$sig=="up",]$SYMBOL
N           <- length(rownames(HCC70_SHELF3_scRNA))#the number of total genes 
K           <- length(cl_up_genes_symbol$SYMBOL)#the number of claudin-low cell lines down-regulated genes  
n           <- length(up_genes)#the number of down-regulated in ELF3-low group 
k           <- length(intersect(cl_up_genes_symbol$SYMBOL, up_genes))  
phyper(q=k-1, m=K, n=N - K, k=n, lower.tail = FALSE)


library(fgsea)
log2fc              <- HCC70_SHELF3_low_VS_high_deres_symbol_highCount
log2fc              <- subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,sig=="up")
rank_vecter         <- log2fc$avg_log2FC
names(rank_vecter)  <- log2fc$SYMBOL
rank_vecter         <- sort(rank_vecter,decreasing = TRUE)
pathway             <- list(cl_up_genes=cl_up_genes_symbol$SYMBOL)
fgseaRes            <- fgsea(pathways = pathway, 
                             stats = rank_vecter,
                             eps = 0.0,scoreType="pos")

# plot the most significantly enriched pathway
plotEnrichment(pathway[[1]],rank_vecter) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway)


####CD24####
####*protein coding genes####
library(readr)
CCLE_log2tpm              <- read_csv("Pro_TNBC/data/CCLE/OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected.csv")
CCLE_log2tpm              <- OmicsExpressionProteinCodingGenesTPMLogp1BatchCorrected
coln                      <- strsplit(colnames(CCLE_log2tpm)[-1],"[ ]")
coln                      <- lapply(coln, function(x) x[1])  %>% unlist()
library(tibble)
CCLE_log2tpm              <- column_to_rownames(CCLE_log2tpm,var = "...1")
CCLE_log2tpm              <- as.matrix(CCLE_log2tpm)  %>% t()
rownames(CCLE_log2tpm)    <- coln
brca_info                 <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
down_gene_id              <- subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,sig=="down")
brca_log2tpm              <- CCLE_log2tpm[,colnames(CCLE_log2tpm) %in% brca_info$DepMap_ID]
brca_log2tpm_down         <- brca_log2tpm[rownames(brca_log2tpm) %in% down_gene_id$SYMBOL,] %>% as.data.frame()
brca_log2tpm_down         <- brca_log2tpm_down %>% as.matrix()  %>% t()
brca_info                 <- arrange(brca_info,lineage_molecular_subtype)
brca_log2tpm_down                   <- brca_log2tpm_down[match(brca_info$DepMap_ID,rownames(brca_log2tpm_down)),]
identical(rownames(brca_log2tpm_down),brca_info$DepMap_ID)
annotation.col              <- brca_info[,22,drop=F] %>% as.data.frame()
rownames(annotation.col)    <- brca_info$DepMap_ID 
colnames(annotation.col)[1] <- "subtype"
annot.colors                <- list(subtype=c("basal_A"="Red","basal_B"= "Blue","HER2_amp"= "Purple","luminal"= "#FFCC00","luminal_HER2_amp"="grey"))
a <- names(table(annotation.col))
library(pheatmap)
p1 <- pheatmap(t(as.matrix(brca_log2tpm_down)),cluster_rows = T,
               cluster_cols = F,
               color = colorRampPalette(c("navy","white","firebrick3"))(100),
               show_colnames = F,border_color = NA,
               scale = "row",
               show_rownames = T,annotation_col   = annotation.col,
               annotation_colors = annot.colors,fontsize=30,cellwidth = 18, cellheight = 25,
               gaps_col =c(13,28,37,48)
)
ggsave(p1,filename = "Pro_TNBC/paper/plot/section_5/heatmap_new.pdf",width = 25,height = 40)


CCLE_log2tpm_df           <- as.data.frame(t(brca_log2tpm))
CCLE_log2tpm_df           <- CCLE_log2tpm_df[match(brca_info$DepMap_ID,rownames(CCLE_log2tpm_df)),]
identical(rownames(CCLE_log2tpm_df),brca_info$DepMap_ID)
CCLE_log2tpm_df$subtype   <- brca_info$lineage_molecular_subtype
ggboxplot(CCLE_log2tpm_df,x = "subtype",y="CD24",add = "point",size= 0.3,add.params=list(size=2))
ggboxplot(CCLE_log2tpm_df,x = "subtype",y="CD44",add = "point",size= 0.3,add.params=list(size=2))


####* all genes####
library(readr)
CCLE_log2tpm              <- read_csv("Pro_TNBC/data/CCLE/OmicsExpressionAllGenesTPMLogp1Profile.csv")
OmicsProfiles             <- subset(OmicsProfiles,Datatype=="rna")
coln                      <- strsplit(colnames(CCLE_log2tpm)[-1],"[ ]")
coln                      <- lapply(coln, function(x) x[1])  %>% unlist()
library(tibble)
CCLE_log2tpm              <- column_to_rownames(CCLE_log2tpm,var = "...1")
CCLE_log2tpm              <- as.matrix(CCLE_log2tpm)  %>% t()
rownames(CCLE_log2tpm)    <- coln
brca_info                 <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
colnames(OmicsProfiles)[3]  <- "DepMap_ID"
brca_info                 <- merge(OmicsProfiles[,c(1,3)],brca_info[,c(3,22)],by="DepMap_ID")
down_gene_id              <- subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,sig=="down")
brca_log2tpm              <- CCLE_log2tpm[,colnames(CCLE_log2tpm) %in%brca_info$ProfileID]
brca_log2tpm_down         <- brca_log2tpm[rownames(brca_log2tpm) %in% down_gene_id$SYMBOL,] %>% as.data.frame()
brca_log2tpm_down         <- brca_log2tpm_down %>% as.matrix()  %>% t()
brca_info                 <- arrange(brca_info,lineage_molecular_subtype)
brca_log2tpm_down                   <- brca_log2tpm_down[match(brca_info$ProfileID,rownames(brca_log2tpm_down)),]
identical(rownames(brca_log2tpm_down),brca_info$ProfileID)
annotation.col              <- brca_info[,3,drop=F] %>% as.data.frame()
rownames(annotation.col)    <- brca_info$ProfileID
colnames(annotation.col)[1] <- "subtype"
annot.colors                <- list(subtype=c("basal_A"="Red","basal_B"= "Blue","HER2_amp"= "Purple","luminal"= "#FFCC00","luminal_HER2_amp"="grey"))
library(pheatmap)
p1 <- pheatmap(t(as.matrix(brca_log2tpm_down)),cluster_rows = T,
               cluster_cols = F,
               color = colorRampPalette(c("navy","white","firebrick3"))(100),
               show_colnames = F,border_color = NA,
               scale = "row",
               show_rownames = T,annotation_col   = annotation.col,
               annotation_colors = annot.colors,fontsize=30,cellwidth = 18, cellheight = 25,
               gaps_col =c(13,28,37,48)
)
ggsave(p1,filename = "Pro_TNBC/paper/plot/section_5/heatmap_new_1.pdf",width = 25,height = 40)


CCLE_log2tpm_df           <- as.data.frame(t(brca_log2tpm))
CCLE_log2tpm_df           <- CCLE_log2tpm_df[match(brca_info$ProfileID,rownames(CCLE_log2tpm_df)),]
identical(rownames(CCLE_log2tpm_df),brca_info$ProfileID)
CCLE_log2tpm_df$subtype   <- brca_info$lineage_molecular_subtype
CCLE_log2tpm_df           <- CCLE_log2tpm_df[,!duplicated(colnames(CCLE_log2tpm_df))]
ggboxplot(CCLE_log2tpm_df,x = "subtype",y="CD24",add = "point",size= 0.3,add.params=list(size=2))
ggboxplot(CCLE_log2tpm_df,x = "subtype",y="CD44",add = "point",size= 0.3,add.params=list(size=2))
save(CCLE_log2tpm_allgenes,file = "Pro_TNBC/data/CCLE/CCLE_log2tpm_allgenes.RData")


Cell_adhesion_molecules    <- c("CLDN4","CDH3","SDC4","VTCN1","MPZL1")
Cell_adhesion_molecules_de <- subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,SYMBOL %in% Cell_adhesion_molecules)



####select Basal-like up genes ####
load("Pro_TNBC/data/CCLE/CCLE.RData")
brca_info                             <- read_csv("Pro_TNBC/output/data/CCLE/brca_info.csv")
brca_info                             <- brca_info[!is.na(brca_info$lineage_molecular_subtype),]
brca_info_BL                          <- subset(brca_info,lineage_molecular_subtype=="basal_A")
brca.log2.rpkm.matrix                 <- CCLE.log2.rpkm.matrix[,colnames(CCLE.log2.rpkm.matrix) %in% brca_info$CCLE_Name]
brca.log2.rpkm.matrix.scale           <- apply(brca.log2.rpkm.matrix,1,scale)
rownames(brca.log2.rpkm.matrix.scale) <- colnames(brca.log2.rpkm.matrix)
brca.log2.rpkm.matrix.BL.scale        <- brca.log2.rpkm.matrix.scale[rownames(brca.log2.rpkm.matrix.scale)%in% brca_info_BL$CCLE_Name,]
brca.log2.rpkm.matrix.BL.scale.median    <- apply(brca.log2.rpkm.matrix.BL.scale,2,median)
brca.log2.rpkm.matrix.BL.scale.median_na <- brca.log2.rpkm.matrix.BL.scale.median[!is.na(brca.log2.rpkm.matrix.BL.scale.median)]
brca.log2.rpkm.matrix.BL.scale.median_na  <- brca.log2.rpkm.matrix.BL.scale.median_na[names(brca.log2.rpkm.matrix.BL.scale.median_na) %in% rownames(HCC70_SHELF3_scRNA)]
plot(sort(brca.log2.rpkm.matrix.BL.scale.median_na))
hist(brca.log2.rpkm.matrix.BL.scale.median_na,breaks = 40)
brca.log2.rpkm.matrix.BL.scale.median_na <- brca.log2.rpkm.matrix.BL.scale.median_na[order(brca.log2.rpkm.matrix.BL.scale.median_na,decreasing = F)]
all(diff(brca.log2.rpkm.matrix.BL.scale.median_na >= 0))
BL_up_genes  <- names(brca.log2.rpkm.matrix.BL.scale.median_na[brca.log2.rpkm.matrix.BL.scale.median_na >0.75])

save(BL_up_genes,file = "Pro_TNBC/paper/data/results/section_5/BL_up_genes.RData")

down_genes  <- HCC70_SHELF3_low_VS_high_deres_symbol_highCount[HCC70_SHELF3_low_VS_high_deres_symbol_highCount$sig=="down",]$ENSEMBL 
N           <- length(rownames(HCC70_SHELF3_scRNA))#the number of total genes 
K           <- length(BL_up_genes)#the number of BLaudin-low cell lines down-regulated genes  
n           <- length(down_genes)#the number of down-regulated in ELF3-low group 
k           <- length(intersect(BL_up_genes, down_genes))  
phyper(q=k-1, m=K, n=N - K, k=n, lower.tail = FALSE)

BL_up_genes_symbol  <- bitr(BL_up_genes,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
library(fgsea)
log2fc              <- HCC70_SHELF3_low_VS_high_deres_symbol_highCount
log2fc              <- subset(HCC70_SHELF3_low_VS_high_deres_symbol_highCount,p_val_adj < 0.05)
rank_vecter         <- log2fc$avg_log2FC
names(rank_vecter)  <- log2fc$SYMBOL
rank_vecter         <- sort(rank_vecter,decreasing = TRUE)
pathway             <- list(BL_up_genes=BL_up_genes_symbol$SYMBOL)
fgseaRes            <- fgsea(pathways = pathway, 
                             stats = rank_vecter,
                             eps = 0.0,scoreType="neg")

# plot the most significantly enriched pathway
plotEnrichment(pathway[[1]],rank_vecter) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway)



brca.log2.rpkm.matrix.BL.scale.median_na.df <- as.data.frame(brca.log2.rpkm.matrix.BL.scale.median_na)
brca.log2.rpkm.matrix.BL.scale.median_na.df  <- brca.log2.rpkm.matrix.BL.scale.median_na.df[order(brca.log2.rpkm.matrix.BL.scale.median_na.df$brca.log2.rpkm.matrix.BL.scale.median_na),,drop=F]
brca.log2.rpkm.matrix.BL.scale.median_na.df$rank  <- rank(brca.log2.rpkm.matrix.BL.scale.median_na.df$brca.log2.rpkm.matrix.BL.scale.median_na)
fig1f_1                              <- ggplot(brca.log2.rpkm.matrix.BL.scale.median_na.df,aes(x=rank,y=brca.log2.rpkm.matrix.BL.scale.median_na))+
  geom_point(size=6,show.legend = F)+ggplot.style+xlab(paste("rank",sep = "\n","(n=22713)"))+ylab("Median expression")+
  scale_x_continuous(limits = c(0, 23000),breaks = seq(0,23000,5000))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.grid.minor = element_blank())+geom_hline(aes(yintercept=0.75,color="red",linewidth=1))
ggsave(fig1f_1,filename = "Pro_TNBC/paper/plot/method/the.rank.charts.of.blgenes.pdf",width=20,height=15)


tumor_cells <- subset(metadata,celltype_major=="Cancer Epithelial")
tumor_scRNA <- scRNA[,colnames(scRNA) %in% tumor_cells$...1]

table(tumor_scRNA$orig.ident)
tumor_cells$Cancer_type <- ifelse(tumor_cells$orig.ident=="CID4523"|tumor_cells$orig.ident=="CID4513","MBC",tumor_cells$subtype)
identical(tumor_cells$...1,colnames(tumor_scRNA))
tumor_scRNA$Cancer_type  <- tumor_cells$Cancer_type
#calculate the proportion of ribosomal genes in cells
tumor_scRNA[["percent.mt"]]     <- PercentageFeatureSet(tumor_scRNA,pattern = "^MT-")
head(tumor_scRNA@meta.data,5)
VlnPlot(tumor_scRNA,features = c("nCount_RNA","nFeature_RNA","percent.mt"),ncol=3)
tumor_scRNA                     <- subset(tumor_scRNA,nFeature_RNA >0 )

#normalized
tumor_scRNA                     <- NormalizeData(tumor_scRNA)

#features choosing
tumor_scRNA                     <- FindVariableFeatures(tumor_scRNA,selection.method = "vst",nfeatures = 2000)#choose high variability between cells.
#identify the 10 most highly variable genes
top10                               <- head(VariableFeatures(tumor_scRNA),10)

#data scaling 
all_gene                               <- rownames(tumor_scRNA)
tumor_scRNA                     <- ScaleData(tumor_scRNA,features = all_gene)#scaled data存放在 pbmc[["RNA"]]@scale.data

#linear dimensionality reduction
tumor_scRNA                     <- RunPCA(tumor_scRNA,features = VariableFeatures(object=tumor_scRNA))

#dimension selection
tumor_scRNA                     <- JackStraw(tumor_scRNA,num.replicate = 100)
tumor_scRNA                     <- ScoreJackStraw(tumor_scRNA,dims = 1:20)
JackStrawPlot(tumor_scRNA,dims=1:20)
ElbowPlot(tumor_scRNA,ndims = 50)#choose dims of 30
#cell clustering
tumor_scRNA                     <- FindNeighbors(tumor_scRNA,dims = 1:40)
tumor_scRNA                     <- FindClusters(tumor_scRNA,resolution = 0.1)
table(Idents(tumor_scRNA))
tumor_scRNA                     <- RunHarmony(tumor_scRNA,"dataset")
#unlinear dimension reduction(UMAP or tSNE)
tumor_scRNA                     <- RunUMAP(tumor_scRNA,dims = 1:40)
tumor_scRNA                     <- RunTSNE(tumor_scRNA,dims = 1:40,check_duplicates = FALSE)
DimPlot(tumor_scRNA,reduction = "umap",label = T)
DimPlot(tumor_scRNA,reduction = "tsne",label = T)


df              <- as.matrix(tumor_scRNA@assays$RNA@data) %>% t() %>% as.data.frame()
df$Cancer_type  <- tumor_scRNA$Cancer_type
ggplot(df,aes(x=Cancer_type,y=KRT5))+geom_boxplot()
ggplot(df,aes(x=Cancer_type,y=KRT14))+geom_boxplot()

FeaturePlot(tumor_scRNA,features = c("KRT5"),#ELF3
            reduction = "umap",pt.size = 1)+ggtitle("KRT5")
FeaturePlot(tumor_scRNA,features = c("KRT14"),#ELF3
            reduction = "umap",pt.size = 1)+ggtitle("KRT14")
Idents(tumor_scRNA)  <- "Cancer_type"
DimPlot(tumor_scRNA,reduction = "umap")
Idents(tumor_scRNA)  <- "orig.ident"
DimPlot(tumor_scRNA,reduction = "umap")
