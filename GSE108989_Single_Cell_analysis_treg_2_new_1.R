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

B_PTR_TTR_Seurat<-CreateSeuratObject(counts = B_PTR_TTR, project = "Treg", min.cells = 3, min.features = 200)
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

library(stringr)
library(Seurat)
library(clusterProfiler)
library(patchwork)
library(dplyr)
library(org.Hs.eg.db)
library(tidyverse)
epimarker<-read.table("D:/GSE139325_T_S_epigenetic_enzyme_for_volcano_plot_new.txt",sep = "\t",header = TRUE,row.names = 1)
for (i in c("wilcox","bimod","roc","t","negbinom","poisson","LR","DESeq2")) {
  if(i == "negbinom" | i == "poisson" | i == "DESeq2"){
    Tregmarker<-FindMarkers(B_PTR_TTR_Seurat,ident.1 = "TTR",ident.2 = "PTR",logfc.threshold = 0,min.pct = 0,test.use = i,slot = "counts")
  }else{
    Tregmarker<-FindMarkers(B_PTR_TTR_Seurat,ident.1 = "TTR",ident.2 = "PTR",logfc.threshold = 0,min.pct = 0,test.use = i)
  }
  
  genename<-rownames(Tregmarker)
  g<-bitr(geneID = genename,fromType = "ENTREZID",toType = "SYMBOL",OrgDb = org.Hs.eg.db,drop = TRUE)
  g<-g[!duplicated(g$SYMBOL),]
  Tregmarker<-Tregmarker[g$ENTREZID,]
  rownames(Tregmarker)<-g$SYMBOL
  Newname<-intersect(toupper(rownames(epimarker)),rownames(Tregmarker))
  Tregmarker_new<-Tregmarker[Newname,]
  write.table(Tregmarker_new,file = paste0("c:/Users/xjmik/Desktop/CRRCmarker",i,".txt"),sep = "\t")
  remove(Tregmarker,Tregmarker_new,Newname,g,genename)
  
  
}
