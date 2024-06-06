if(T){
  library(patchwork)
  library(Seurat)#4.0
  library(tidyverse)
  library(SeuratData)
  library(cowplot)
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(clustree)
  library(scDblFinder)
  library(BiocParallel)
  library(SingleCellExperiment)
  library(scater)
  library(CellChat)
  library(ggplot2)
  library(ggalluvial)
  library(VGAM)
  library(nichenetr)
  library(GSVA)
  library(qusage)
  library(Hmisc)
  library(pheatmap)
  library(SCENIC)
  library(hdf5r)
  library(MAST)
  library(devtools)
  library(harmony)
  library(future)
  library(glmGamPoi)
  library(COSG)
  library(job)
  library(reshape2)
  library(ggpubr)
}
setwd("/sepsis/")
expr=read.csv("scp_gex_matrix_raw.csv")
head(expr[, 1:5], n = 5)
rownames(expr)<-expr[,1]
expr <- expr[,-1]
expr <- as.matrix(expr)
scRNA <- CreateSeuratObject(counts = expr,min.cells = 1)
scRNA$barcode=rownames(scRNA@meta.data)
list=read.table("scp_meta_change.txt",header = T)
scRNA@meta.data=left_join(scRNA@meta.data, list, by='barcode')
rownames(scRNA@meta.data)=scRNA@meta.data$barcode

####scDblFinder####
scRNA = as.SingleCellExperiment(scRNA)
scRNA <- scDblFinder(scRNA, samples="orig.ident", BPPARAM=MulticoreParam(8))
scRNA=as.Seurat(scRNA)
####QC####
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
pdf("p1.pdf")
scRNA@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nCount_RNA, fill=orig.ident)) + 
  geom_density(alpha = 0.15) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("log10 Counts density") +
  geom_vline(xintercept = 1000) +
  geom_vline(xintercept = 10000)
dev.off()
pdf("p2.pdf")
scRNA@meta.data %>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() +
  theme_classic() +
  ylab("log10 gene density") +
  geom_vline(xintercept = 500)+
  geom_vline(xintercept = 600)+
  geom_vline(xintercept = 2000)+
  geom_vline(xintercept = 3000)
dev.off()
pdf("p3.pdf")
scRNA@meta.data %>% 
  ggplot(aes(color=orig.ident, x=percent.mt, fill=orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  geom_vline(xintercept = 5) +
  geom_vline(xintercept = 10)
dev.off()

scRNA = subset(scRNA,cells=which(scRNA$nFeature_RNA > 500 & scRNA$nFeature_RNA < 3000 & scRNA$percent.mt < 5 
& scRNA$nCount_RNA > 1000 & scRNA$nCount_RNA < 10000 & scRNA$scDblFinder.class=="singlet"))
####Normalize Data####
scRNA <- NormalizeData(scRNA, normalization.method = 'LogNormalize', scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA, vars.to.regress = "percent.mt", features = VariableFeatures(scRNA))
scRNA <- RunPCA(scRNA)
####Removing batch effects####
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony=20,reduction = "pca");gc()
pct<-scRNA[["harmony"]]@stdev/sum(scRNA[["harmony"]]@stdev)*100
cumu<-cumsum(pct)
co1<-which(cumu >80 & pct<5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 <- min(co1,co2)
co2 
pc.num=1:co2
#RunUMAP,RunTSNE
scRNA <- RunUMAP(scRNA,reduction="harmony", dims=pc.num,seed.use=123456L)
scRNA <- RunTSNE(scRNA,reduction="harmony", dims=pc.num,seed.use=123456L)
scRNA <- FindNeighbors(scRNA,reduction="harmony",dims = pc.num,k.param = 20)
scRNA <- FindClusters(scRNA, resolution = seq(0.5,1,by=0.1))
ggsave('clustree.pdf', clustree(scRNA), width = 12, height = 10)
scRNA$seurat_clusters=scRNA$RNA_snn_res.0.8
saveRDS(scRNA, "Cluster_0.8.RDS")
Idents(scRNA)=scRNA$seurat_clusters
ggsave('DimPlot.pdf', DimPlot(scRNA,label = T), width = 12, height = 10)
####Cell annotation####
#first round
p=VlnPlot(scRNA,features = c("CD3D","CD3E","CD3G","CD14","CD4","CD8A","CD79A","CLEC4C","NCAM1","NKG7","GNLY","KLRD1","ITGA2B","GP9","CEACAM8"),pt.size = 0,stack = T)
ggsave('VlnPlot.big.pdf', p, width = 12, height = 10)
scRNA = RenameIdents(scRNA,
                     "0"="Mono/Macrophages",
                     "1"="T cell",
                     "2"="Mono/Macrophages",
                     "3"="T cell",
                     "4"="Mono/Macrophages",
                     "5"="NK cell",
                     "6"="pDC",
                     "7"="DC",
                     "8"="NK cell",
                     "9"="Mono/Macrophages",
                     "10"="B cell",
                     "11"="B cell",
                     "12"="Mono/Macrophages",
                     "16"="Mono/Macrophages"
)
scRNA$Cell_type=Idents(scRNA)
saveRDS(scRNA, "Cluster_0.8_Celltype.RDS")
p=DimPlot(scRNA,label = T,reduction = 'tsne')
ggsave('1.cellType.pdf', p, width = 7, height = 5)
#chart
markers = c("CD14","VCAN","FCN1",#Myeloid
  "CD3D","CD3E","CD3G",#T
  "NCAM1","NKG7","GNLY","KLRD1",#NK
  "CLEC4C","JCHAIN","IL3RA","LILRA4",#pDC
  "CD1C","FCER1A","CLEC10A",#DC
 "CD79A","CD19","MS4A1")#B
p <- DotPlot(object = scRNA, features=markers, group.by= "Cell_type",cols = c("lightgrey", "red"))+
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
     geom_vline(xintercept = c(3.5,6.5,10.5,14.5,17.5), linetype = "dashed", color = "lightgrey")+
     geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5,5.5), linetype = "dashed", color = "lightgrey") 
ggsave('2.markerDotPlot.pdf', p, width = 8, height = 4)
#################healthy，disease donors#################
#Cell ratio bar diagram
scRNA@meta.data$donor <- factor(scRNA@meta.data$donor, levels = c("P18F","P17H","P20H","P15F","P08H","P13H","P07H","P06F","P04H","C2P01H","P09H","P02H","C2P05F","C2P07H","C2P13F","C2P16H","C2P10H","C2P19H","C2P15H","P08","P04","P09","C2P02","C2P01","C2P09","C2P05","P05","C2P13","C2P11","C2P12","C2P10","C2P15","C2P16","C2P19","C2P18","C2P21","P633","P636","P640","P662","P669","P670","P672","P671","E25","E1","E12","E16"))
p1=scRNA@meta.data %>%
  ggplot(aes(x=factor(donor),fill=Cell_type))+
  geom_bar(position = position_fill())+
  scale_fill_manual(values=c("#A2C8D9","#ACD486","#36973B",
                                      "#E49B9B",
                                      "#F6BE73","#CCB5D3","#683E9C",
                                      "#FED703","#B25A28",
                                      "#D9D9D9",
                                      "#FCE1D4"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("")+ylab(paste0("Proportion of cells"))
#cell number
p2=scRNA@meta.data %>%
  ggplot(aes(x=disease,fill=Cell_type))+
  geom_bar(position = position_fill())+
  scale_fill_manual(values=c("#A2C8D9","#ACD486","#36973B",
                                      "#E49B9B",
                                      "#F6BE73","#CCB5D3","#683E9C",
                                      "#FED703","#B25A28",
                                      "#D9D9D9",
                                      "#FCE1D4"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("")+ylab(paste0("Proportion of cells"))
combined_plot <- plot_grid(p1, p2, ncol = 2, rel_widths = c(0.6, 0.4))
ggsave('3.cell-proportion.pdf', plot = combined_plot, width = 10, height = 4)
#cell ratio
Cellratio <- as.data.frame(prop.table(table(scRNA@meta.data$Cell_type, scRNA$donor), margin = 2))
colnames(Cellratio)= c("cell","donor","Freq")
list=scRNA@meta.data[,c(6,7)]
Cellratio=left_join(Cellratio, list, by='donor')
data=Cellratio[,c(4,1,3)]
colnames(data)=c("Group","Cell","Freq")
p=ggboxplot(data, x="Cell", y="Freq", color = "Group", 
            ylab="",xlab="",
            legend.title="Group",
            palette = c("#9cbfdb","#8bad83"),
            width=0.6, add = "none")
p1=p+stat_compare_means(aes(group=Group),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                        label = "p.signif")
ggsave('3.cell-proportion-difference.pdf', p1, width = 7, height = 4)
#SMR result + scRNA
gene=c('LINC00211','SPATC1','REV1','C4B','C4A','ABCC1','ARL14EP','SNORD89','MANEA','LGALS9','ORM1','ITIH2','SLC35E2B','SIGLEC7','PSPH','ORM2','MIB2')
  gene_matrix=scRNA@assays$RNA@data[which(row.names(scRNA@assays$RNA@data)%in%gene),]
  gene_matrix=t(gene_matrix)
  gene_matrix=as.data.frame(gene_matrix)
  gene_matrix$barcode=row.names(gene_matrix)
  list=scRNA@meta.data[,c(4,6)]
  gene_matrix=left_join(gene_matrix, list, by='barcode')
  rownames(gene_matrix)=gene_matrix$barcode
  gene_matrix=gene_matrix[,c(15,1:13)]
  #chart
  colnames(gene_matrix)[1]="Group"
  data=melt(gene_matrix,id.vars=c("Group"))
  colnames(data)=c("Group","Gene","Expression")
  data <- data[data$Expression > 0, ]
  p=ggboxplot(data, x="Gene", y="Expression", color = "Group", 
              ylab=paste0("Genes in ",cell[i]),xlab="",
              legend.title="Group",
              palette = c("#9cbfdb","#8bad83"),
              width=0.6, add = "none")
  p1=p+stat_compare_means(aes(group=Group),
                          method="wilcox.test",
                          symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                          label = "p.signif")+
    theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust=1))
  ggsave(paste0("5.gene-difference-allcell.pdf"),p1, width = 6, height = 4)
#################T cell subtype#################
cell_type = "T cell"
aim_cell = subset(scRNA,idents=cell_type)
aim_cell <- NormalizeData(aim_cell, normalization.method = 'LogNormalize', scale.factor = 10000)
aim_cell <- FindVariableFeatures(aim_cell, selection.method = "vst", nfeatures = 2000)
aim_cell <- ScaleData(aim_cell, vars.to.regress = "percent.mt", features = VariableFeatures(aim_cell))
aim_cell <- RunPCA(aim_cell)
aim_cell <- RunHarmony(aim_cell, group.by.vars="donor", max.iter.harmony=20,reduction = "pca");gc()
ElbowPlot(aim_cell)
pct<-aim_cell[["harmony"]]@stdev/sum(aim_cell[["harmony"]]@stdev)*100
cumu<-cumsum(pct)
co1<-which(cumu >80 & pct<5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 <- min(co1,co2)
co2
pc.num=1:co2
aim_cell <- RunUMAP(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- RunTSNE(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- FindNeighbors(aim_cell,reduction="harmony",dims = pc.num,k.param = 20)
aim_cell <- FindClusters(aim_cell,resolution = seq(0.1,0.5,by=0.1))
clustree(aim_cell)
aim_cell$seurat_clusters=aim_cell$RNA_snn_res.0.3
Idents(aim_cell)=aim_cell$seurat_clusters
DimPlot(aim_cell,label = T,reduction = "tsne")
#marker
markers = c("CD4","CD8A","CD8B",
            "MKI67","TOP2A"#Cycling T cells
)
DotPlot(object = aim_cell, features=markers, group.by= "seurat_clusters",cols = c("lightgrey", "red"))
aim_cell = RenameIdents(aim_cell,
                        "0"="CD8 T cell",
                        "1"="CD4 T cell",
                        "2"="CD4 T cell",
                        "3"="CD4 T cell",
                        "4"="CD8 T cell"
)
aim_cell$subType=Idents(aim_cell)
DimPlot(aim_cell,label = T,reduction = "tsne")
p=DimPlot(aim_cell,label = T,reduction = "tsne")
ggsave("7.1.T-celltype-TSNE.pdf",p, width = 9, height = 6)
p=DotPlot(object = aim_cell, features=markers, group.by= "subType",cols = c("lightgrey", "red"))+
     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("7.2.T-celltype-DotPlot.pdf",p, width = 6, height = 4)
saveRDS(aim_cell, "Cluster_0.8_subtype_T.RDS")
########CD4 T cell subtype########
cell_type = "CD4 T cell"
aim_cell2 = subset(aim_cell,idents=cell_type)
aim_cell2 <- NormalizeData(aim_cell2, normalization.method = 'LogNormalize', scale.factor = 10000)
aim_cell2 <- FindVariableFeatures(aim_cell2, selection.method = "vst", nfeatures = 2000)
aim_cell2 <- ScaleData(aim_cell2, vars.to.regress = "percent.mt", features = VariableFeatures(aim_cell2))
aim_cell2 <- RunPCA(aim_cell2)
aim_cell2 <- RunHarmony(aim_cell2, group.by.vars="donor", max.iter.harmony=20,reduction = "pca");gc()
ElbowPlot(aim_cell2)
pct<-aim_cell2[["harmony"]]@stdev/sum(aim_cell2[["harmony"]]@stdev)*100
cumu<-cumsum(pct)
co1<-which(cumu >80 & pct<5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 <- min(co1,co2)
co2
pc.num=1:co2
aim_cell2 <- RunUMAP(aim_cell2,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell2 <- RunTSNE(aim_cell2,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell2 <- FindNeighbors(aim_cell2,reduction="harmony",dims = pc.num,k.param = 20)
aim_cell2 <- FindClusters(aim_cell2,resolution = seq(0.1,0.5,by=0.1))
clustree(aim_cell2)
aim_cell2$seurat_clusters=aim_cell2$RNA_snn_res.0.4
Idents(aim_cell2)=aim_cell2$seurat_clusters
DimPlot(aim_cell2,label = T,reduction = "tsne")
#all marker
markers = c(
            "CCR7","SELL","TCF7","LEF1","TXK",#Tn,Naive T cell
            "CXCR5",# CXCR5+ pre-Tfh,Follicular helper T cell
            "ADSL",#ADSL+ Tn
            "IL7R",#IL7R- Tn,IL7R是memory的marker
            "TNF",#TNF+ T
            "AGER",#AGER+ Tm,memory T cell
            "TIMP1",#TIMP1+ Tm
            "CREM",#CREM+ Tm,CPAG+CREM- Tm
            "CCL5",#CCL5+ Tm
            #"CPAG",#CPAG+ Tm,CCL5+
            "GZMK","KLRG1",#GZMK+ Tem,Effector memory T Cell,CCL5+CPAG-
            "CX3CR1","TBX21",#Temra
            "RORA","RORC",#Th17
            "CCR6",#CCR6+ Th17
            "IL26",#IL26+ Th17
            "TOX","TOX2","IL21","CXCL13","GNG4",#IL21+ Tfh
            "IFNG",#IFNG+ Tfh
            "RTKN2","IL2RA",#Treg
            "TNFRSF9",#TNFRSF9- Treg,S1PR1-;TNFRSF9+ Treg,S1PR1-
            "S1PR1",#S1PR1+ Treg,TNFRSF9-
            "STAT1","IFIT1","IRF7",#ISG+
            "NME1","NME2","CCR4" #NME1+ CCR4- T,NME1+ CCR4+ T
)
markers = c(
  "CCR7","TCF7",#Tn,Naive T cell
  "IL7R","TIMP1",#TIMP1+ Tm
  "RTKN2","IL2RA","STAT1",#ISG+ Treg
  "CCL5","GZMK","KLRG1" #GZMK+ Tem
)

DotPlot(object = aim_cell2, features=markers, group.by= "seurat_clusters",cols = c("lightgrey", "red"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
aim_cell2 = RenameIdents(aim_cell2,
                         "0"="Tn", 
                         "1"="TIMP1+ Tm", #IL7R,TIMP1
                         "2"="Tn", #CCR7,SELL,TCF7,LEF1,TXK
                         "3"="Tn", 
                         "4"="ISG+ Treg",#RORA,RORC,CCR6
                         "5"="GZMK+ Tem" #CCL5,GZMK,KLRG1
)
aim_cell2$subType=Idents(aim_cell2)
DimPlot(aim_cell2,label = T,reduction = "tsne")
p=DimPlot(aim_cell2,label = T,reduction = "tsne")
ggsave("10.1.CD4T-celltype-TSNE.pdf",p, width = 6, height = 4)
p=DotPlot(object = aim_cell2, features=markers, group.by= "subType",cols = c("lightgrey", "red"))+
  theme(axis.text.x = element_text(angle = 35, vjust = 0.5, hjust=1))+
  geom_vline(xintercept = c(2.5,4.5,7.5), linetype = "dashed", color = "lightgrey")+
  geom_hline(yintercept = c(1.5, 2.5,3.5), linetype = "dashed", color = "lightgrey")
ggsave("10.2.CD4T-celltype-DotPlot.pdf",p, width = 8, height = 4)
#CD4T cell subtype radio
Cellratio <- as.data.frame(prop.table(table(aim_cell2@meta.data$subType, aim_cell2$donor), margin = 2))
colnames(Cellratio)= c("cell","donor","Freq")
list=aim_cell2@meta.data[,c(6,7)]
Cellratio=left_join(Cellratio, list, by='donor')
data=Cellratio[,c(4,1,3)]
colnames(data)=c("Group","Cell","Freq")
p=ggboxplot(data, x="Cell", y="Freq", color = "Group", 
            ylab="",xlab="",
            legend.title="Group",
            palette = c("#9cbfdb","#8bad83"),
            width=0.6, add = "none")
p1=p+stat_compare_means(aes(group=Group),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                        label = "p.signif")

ggsave("10.3.CD4T-cell-proportion-wilcox.test.pdf",p1, width = 6, height = 4)
saveRDS(aim_cell2, "Cluster_0.8_subtype_CD4T.RDS")
p=aim_cell2@meta.data %>%
  ggplot(aes(x=disease,fill=subType))+
  geom_bar(position = position_fill())+
  scale_fill_manual(values=c("#A2C8D9","#ACD486","#36973B",
                             "#E49B9B",
                             "#F6BE73","#CCB5D3","#683E9C",
                             "#FED703","#B25A28",
                             "#D9D9D9",
                             "#FCE1D4"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("")+ylab(paste0("Proportion of cells"))
ggsave("10.3.CD4T-cell-proportion.pdf",p, width = 4, height = 4)

p1 <- FeaturePlot(object = aim_cell2, features = "LGALS9",reduction = "tsne",cols = c("lightgrey", "#7301A8FF"))
ggsave('10.4.CD4T-FeaturePlot-LGALS9.pdf', p1, width = 8, height = 6)

#################Myeloid cell subtype#################
cell_type = "Mono/Macrophages"
aim_cell = subset(scRNA,idents=cell_type)
aim_cell <- NormalizeData(aim_cell, normalization.method = 'LogNormalize', scale.factor = 10000)
aim_cell <- FindVariableFeatures(aim_cell, selection.method = "vst", nfeatures = 2000)
aim_cell <- ScaleData(aim_cell, vars.to.regress = "percent.mt", features = VariableFeatures(aim_cell))
aim_cell <- RunPCA(aim_cell)
aim_cell <- RunHarmony(aim_cell, group.by.vars="donor", max.iter.harmony=20,reduction = "pca");gc()
ElbowPlot(aim_cell)
pct<-aim_cell[["harmony"]]@stdev/sum(aim_cell[["harmony"]]@stdev)*100
cumu<-cumsum(pct)
co1<-which(cumu >80 & pct<5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
co2 <- min(co1,co2)
co2
pc.num=1:co2
#降维聚类
aim_cell <- RunUMAP(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- RunTSNE(aim_cell,reduction="harmony", dims=pc.num,seed.use=123456L)
aim_cell <- FindNeighbors(aim_cell,reduction="harmony",dims = pc.num,k.param = 20)
aim_cell <- FindClusters(aim_cell,resolution = seq(0.1,0.5,by=0.1))
clustree(aim_cell)
aim_cell$seurat_clusters=aim_cell$RNA_snn_res.0.3
Idents(aim_cell)=aim_cell$seurat_clusters
DimPlot(aim_cell,label = T,reduction = "tsne")
#all marker
markers = c("FAM129C","GZMB","IGJ","PTPRCAP","LILRA4","IGKC","IRF4",#pDC_LILRA4
            #"XCR1","CLEC9A","CADM1","SLAMF8",#cDC1_CLEC9A
            "CD1C","CLEC10A","CD1E",#cDC2_CD1C
            #"LAMP3","CCR7","CST7","IL4I1",#cDC3_LAMP3
            "S100A8","S100A9","FCN1",#Mono_CD14
            "HK3","FCGR3A","PILRA",#Mono_CD16
            "VCAN","TREM1","OLR1","VEGFA",#Macro_VCAN
            "CXCL10","ISG15","CCL2","ILIRN",#Macro_ISG15
            "APOE","RNASE1","C1QC","APOC1","SLCO2B1","GPNMB","TREM2"#Macro_C1QC
)
markers = c(
            "S100A8","S100A9","FCN1",#Mono_CD14
            "HK3","FCGR3A","PILRA",#Mono_CD16
            "ISG15",#Macro_ISG15
            "C1QC"#Macro_C1QC
)
DotPlot(object = aim_cell, features=markers, group.by= "seurat_clusters",cols = c("lightgrey", "red")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
aim_cell = RenameIdents(aim_cell,
                        "0"="Mono_CD14",
                        "1"="Mono_CD14",
                        "2"="Mono_CD16",
                        "3"="Macro_ISG15",
                        "4"="Mono_CD14",
                        "5"="Macro_C1QC"
)
aim_cell$subType=Idents(aim_cell)
DimPlot(aim_cell,label = T,reduction = "tsne")
p=DimPlot(aim_cell,label = T,reduction = "tsne")
ggsave("8.1.Myeloid-celltype-TSNE.pdf",p, width = 7, height = 4)
p=DotPlot(object = aim_cell, features=markers, group.by= "subType",cols = c("lightgrey", "red"))+
     theme(axis.text.x = element_text(angle = 25, vjust = 0.5, hjust=1))+
  geom_vline(xintercept = c(3.5,6.5,7.5), linetype = "dashed", color = "lightgrey")+
  geom_hline(yintercept = c(1.5, 2.5,3.5), linetype = "dashed", color = "lightgrey") 
ggsave("8.2.Myeloid-celltype-DotPlot.pdf",p, width = 6, height = 4)
saveRDS(aim_cell, "Cluster_0.8_subtype_Myeloid.RDS")

Cellratio <- as.data.frame(prop.table(table(aim_cell@meta.data$subType, aim_cell$donor), margin = 2))
colnames(Cellratio)= c("cell","donor","Freq")
list=aim_cell@meta.data[,c(6,7)]
Cellratio=left_join(Cellratio, list, by='donor')
data=Cellratio[,c(4,1,3)]
colnames(data)=c("Group","Cell","Freq")
p=ggboxplot(data, x="Cell", y="Freq", color = "Group", 
            ylab="",xlab="",
            legend.title="Group",
            palette = c("#9cbfdb","#8bad83"),
            width=0.6, add = "none")

p1=p+stat_compare_means(aes(group=Group),
                        method="wilcox.test",
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                        label = "p.signif")
ggsave("8.3.Myeloid-cell-proportion-wilcox.test.pdf",p1, width = 6, height = 4)

aim_cell@meta.data$donor <- factor(aim_cell@meta.data$donor, levels = c("P18F","P17H","P20H","P15F","P08H","P13H","P07H","P06F","P04H","C2P01H","P09H","P02H","C2P05F","C2P07H","C2P13F","C2P16H","C2P10H","C2P19H","C2P15H","P08","P04","P09","C2P02","C2P01","C2P09","C2P05","P05","C2P13","C2P11","C2P12","C2P10","C2P15","C2P16","C2P19","C2P18","C2P21","P633","P636","P640","P662","P669","P670","P672","P671","E25","E1","E12","E16"))
p1=aim_cell@meta.data %>%
  ggplot(aes(x=factor(donor),fill=subType))+
  geom_bar(position = position_fill())+
  scale_fill_manual(values=c("#A2C8D9","#ACD486","#36973B",
                             "#E49B9B",
                             "#F6BE73","#CCB5D3","#683E9C",
                             "#FED703","#B25A28",
                             "#D9D9D9",
                             "#FCE1D4"))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("")+ylab(paste0("Proportion of cells"))
p2=aim_cell@meta.data %>%
  ggplot(aes(x=disease,fill=subType))+
  geom_bar(position = position_fill())+
  scale_fill_manual(values=c("#A2C8D9","#ACD486","#36973B",
                                      "#E49B9B",
                                      "#F6BE73","#CCB5D3","#683E9C",
                                      "#FED703","#B25A28",
                                      "#D9D9D9",
                                      "#FCE1D4"))+
  theme_classic() +
  xlab("")+ylab(paste0("Proportion of cells"))
combined_plot <- plot_grid(p1, p2, ncol = 2, rel_widths = c(0.6, 0.4))
ggsave("8.3.Myeloid-cell-proportion.pdf",combined_plot, width = 8, height = 4)

p1 <- FeaturePlot(object = aim_cell, features = "LGALS9",reduction = "tsne",cols = c("lightgrey", "#7301A8FF"))
ggsave('8.4.Myeloid-FeaturePlot-LGALS9.pdf', p1, width = 8, height = 6)

gene=c("LGALS9")
p1=DotPlot(object = scRNA, features=gene, group.by= "Cell_type",cols = c("lightgrey", "red"))
p2=DotPlot(object = aim_cell, features=gene, group.by= "subType",cols = c("lightgrey", "red"))
ggsave("8.4.Myeloid-DotPlot-LGALS9.pdf",p, width = 5, height = 4)

##############pseudotime analysis##############
library(monocle3)
aim_cells=readRDS("Cluster_0.8_subtype_CD4T.RDS")
fData <- data.frame(gene_short_name = row.names(aim_cell@assays$RNA@counts),
                    row.names = row.names(aim_cell@assays$RNA@counts))
cds <- new_cell_data_set(aim_cell@assays$RNA@counts,
                         cell_metadata = aim_cell@meta.data,
                         gene_metadata = fData)
cds <- preprocess_cds(cds)
plot_pc_variance_explained(cds)
cds <- preprocess_cds(cds, num_dim =10,verbose=T)
cds <- reduce_dimension(cds,reduction_method  = "UMAP")

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(aim_cell, reduction = "tsne")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed

plot_cells(cds,color_cells_by="subType")
cds <- cluster_cells(cds,cluster_method = "louvain")
cds <- learn_graph(cds)
cds <- order_cells(cds)
p=plot_cells(cds,color_cells_by="subType",label_branch_points=F,label_leaves=F)+
  plot_cells(cds,color_cells_by="pseudotime",label_branch_points=F,label_leaves=F)
ggsave("10.7.CD4T-Tn-pseudotime.pdf",p, width = 10, height = 4)

#################cellchat#################
sc <- readRDS("Cluster_0.8_subtype_Myeloid.RDS")
cellchat <- createCellChat(object = sc,
                           meta = sc@meta.data,
                           group.by = "subType")
CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# set the used database in the object
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 8) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# project gene expression data onto PPI network (optional)
# cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

#all the inferred cell-cell communications at the level of ligands/receptors
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "cell-cell_communications.all.csv")

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
p=netVisual_bubble(cellchat, sources.use = c(1),targets.use = c(1,2,3,4), remove.isolate = FALSE)
ggsave("10.8.CD4T.Tn.Cellchat.pdf",p,width = 6, height = 4)

#############GO、KEGG#############
library(org.Hs.eg.db)
list=read.table("eqtlgen-gene.txt",header = T)
entrezid_all = mapIds(x = org.Hs.eg.db,
                      keys = list$gene,
                      keytype = "SYMBOL",
                      column = "ENTREZID")
entrezid_all  = na.omit(entrezid_all)
entrezid_all = data.frame(entrezid_all)
head(entrezid_all)
###GO###
GO_enrich = enrichGO(gene = entrezid_all[,1],
                     OrgDb = org.Hs.eg.db, 
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pvalueCutoff = 0.05,qvalueCutoff = 0.05,
                     readable = T)

GO_enrich  = data.frame(GO_enrich) 
write.csv(GO_enrich,'GO_enrich.csv')

library(readr)
go_enrich=read_csv("GO_enrich.csv")
x=go_enrich$Count
y=factor(go_enrich$Description,levels = go_enrich$Description)#设置Y轴
p = ggplot(go_enrich,aes(x,y))+
  geom_point() +
  theme_bw() + 
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_line(color = "transparent"),
        panel.grid.minor = element_line(color = "transparent"))
p1 = p + geom_point(aes(size=Count,color=-0.5*log(p.adjust),shape=ONTOLOGY,))+
  scale_color_gradient(low = "gold", high = "red")
p2 = p1 + labs(color=expression(-log[10](p.adjust)),
               size="Count",
               x="Count",
               y="Go_term",
               title="Go enrichment of test Genes")+scale_x_continuous(breaks=seq(min(go_enrich$Count), max(go_enrich$Count), by=1))
ggsave('GO.pdf', p2, width = 14, height = 7)

###KEGG###
KEGG_enrich = enrichKEGG(gene = entrezid_all[,1],
                         keyType = "kegg",
                         pAdjustMethod = 'fdr',
                         organism= "hsa",
                         qvalueCutoff = 0.05,
                         pvalueCutoff=0.05)
KEGG_enrich  = data.frame(KEGG_enrich)
write.csv(KEGG_enrich,'KEGG_enrich.csv')

