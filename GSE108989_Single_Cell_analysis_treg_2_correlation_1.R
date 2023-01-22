A<-read.table("D:/Jmjd1c_Treg_Tumor/GSE108989/GSE108989_CRC.TCell.S10805.norm.centered.txt",header = TRUE,sep = "\t")
B<-read.table("D:/Jmjd1c_Treg_Tumor/GSE108989/GSE108989_CRC.TCell.S11138.count.txt",header = TRUE,sep = "\t")
C<-read.table("D:/Jmjd1c_Treg_Tumor/GSE108989/GSE108989_CRC.TCell.S11138.TPM.txt",header = TRUE,sep = "\t")

library(dplyr)

A_PTR_TTR<-select(A,starts_with("geneID"),starts_with("geneSymbol"),starts_with("PTR"),starts_with("TTR"))

B_PTR_TTR<-select(B,starts_with("geneID"),starts_with("symbol"),starts_with("PTR"),starts_with("TTR"))

C_PTR_TTR<-select(C,starts_with("geneID"),starts_with("symbol"),starts_with("PTR"),starts_with("TTR"))

A_PTR<-select(A_PTR_TTR,starts_with("PTR"))
A_TTR<-select(A_PTR_TTR,starts_with("TTR"))

B_PTR<-select(B_PTR_TTR,starts_with("PTR"))
B_TTR<-select(B_PTR_TTR,starts_with("TTR"))

C_PTR<-select(C_PTR_TTR,starts_with("PTR"))
C_TTR<-select(C_PTR_TTR,starts_with("TTR"))

rowname_A<-A[,1]
rowname_B<-B[,1]
rowname_C<-C[,1]

rownames(A_PTR)<-rowname_A
rownames(A_TTR)<-rowname_A
rownames(B_PTR)<-rowname_B
rownames(B_TTR)<-rowname_B
rownames(C_PTR)<-rowname_C
rownames(C_TTR)<-rowname_C

library(cowplot)
library(Seurat)


rownames(B_PTR_TTR)<-rowname_B
B_PTR_TTR<-B_PTR_TTR[,-1]
B_PTR_TTR<-B_PTR_TTR[,-1]
metadata<-read.table("D:/Jmjd1c_Treg_Tumor/GSE108989/GSE108989_metadata.txt",sep = "\t",header = TRUE)
PTR_metadata<-filter(metadata,sampleType == "PTR")
TTR_metadata<-filter(metadata,sampleType == "TTR")
PTR_TTR_metadata<-rbind(PTR_metadata,TTR_metadata)
Z<-as.data.frame(colnames(B_PTR_TTR))
rownames(PTR_TTR_metadata)<-PTR_TTR_metadata[,1]
B_PTR_TTR_Seurat<-CreateSeuratObject(counts = B_PTR_TTR, project = "Treg", min.cells = 3, min.features = 200,meta.data = PTR_TTR_metadata)
B_PTR_TTR_Seurat <- subset(B_PTR_TTR_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
B_PTR_TTR_Seurat <- NormalizeData(B_PTR_TTR_Seurat)
B_PTR_TTR_Seurat <- FindVariableFeatures(B_PTR_TTR_Seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(B_PTR_TTR_Seurat)
B_PTR_TTR_Seurat <- ScaleData(B_PTR_TTR_Seurat, features = all.genes)
B_PTR_TTR_Seurat <- RunPCA(B_PTR_TTR_Seurat, features = VariableFeatures(object = B_PTR_TTR_Seurat))
ElbowPlot(B_PTR_TTR_Seurat)
B_PTR_TTR_Seurat <- FindNeighbors(B_PTR_TTR_Seurat, dims = 1:10)
B_PTR_TTR_Seurat <- FindClusters(B_PTR_TTR_Seurat, resolution = 0.5)

PTR<-colnames(B_PTR)
TTR<-colnames(B_TTR)

Idents(B_PTR_TTR_Seurat, cells = PTR)<-'PTR'
Idents(B_PTR_TTR_Seurat, cells = TTR)<-'TTR'

B_PTR_TTR_Seurat <- RunUMAP(B_PTR_TTR_Seurat, dims = 1:10)
B_PTR_TTR_Seurat <- RunTSNE(B_PTR_TTR_Seurat, dims = 1:10)

DimPlot(B_PTR_TTR_Seurat, reduction = "umap")
DimPlot(B_PTR_TTR_Seurat, reduction = "tsne")

B_TTR_seurat<-subset(B_PTR_TTR_Seurat,idents = "TTR")

Idents(B_TTR_seurat)<-B_TTR_seurat@meta.data$Patient_ID

for (i in unique(B_TTR_seurat@meta.data[["Patient_ID"]])) {
  assign(i,subset(B_TTR_seurat,idents = i))
  assign("a",slot(get(i),"assays"))
  assign("b",as.data.frame(t(as.data.frame(a$RNA@scale.data))))
  assign("c",select(b,one_of("221037","3458")))
  colnames(c)<-c("JMJD1C","IFNG")
  assign(paste("scaldata_",i,"_TTR",sep = ""),data.frame(mean(c$JMJD1C),mean(c$IFNG)))
}

scaldata<-rbind(scaldata_P0123_TTR,scaldata_P0215_TTR,scaldata_P0309_TTR,scaldata_P0411_TTR,scaldata_P0413_TTR,scaldata_P0825_TTR,scaldata_P0909_TTR,scaldata_P1012_TTR,scaldata_P1212_TTR,scaldata_P1228_TTR)
colnames(scaldata)<-c("JMJD1C","IFNG")
library(ggplot2)
library(ggpubr)
ggplot(data = scaldata,aes(x=JMJD1C,y=IFNG)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=TRUE,size=1.5,color="red")+stat_cor(data = scaldata,method = "spearman") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
