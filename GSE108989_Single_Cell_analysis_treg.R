A<-read.table("G:/GSE108989/GSE108989_CRC.TCell.S10805.norm.centered.txt",header = TRUE,sep = "\t")
B<-read.table("G:/GSE108989/GSE108989_CRC.TCell.S11138.count.txt",header = TRUE,sep = "\t")
C<-read.table("G:/GSE108989/GSE108989_CRC.TCell.S11138.TPM.txt",header = TRUE,sep = "\t")

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

B_PTR_seurat <- CreateSeuratObject(counts = B_PTR, project = "IMMUNE_CTRL", min.cells = 5)
B_PTR_seurat$stim <- "CTRL"
B_PTR_seurat <- subset(B_PTR_seurat, subset = nFeature_RNA > 500)
B_PTR_seurat <- NormalizeData(B_PTR_seurat, verbose = FALSE)
B_PTR_seurat <- FindVariableFeatures(B_PTR_seurat, selection.method = "vst", nfeatures = 2000)

B_TTR_seurat <- CreateSeuratObject(counts = B_TTR, project = "IMMUNE_STIM", min.cells = 5)
B_TTR_seurat$stim <- "STIM"
B_TTR_seurat <- subset(B_TTR_seurat, subset = nFeature_RNA > 500)
B_TTR_seurat <- NormalizeData(B_TTR_seurat, verbose = FALSE)
B_TTR_seurat <- FindVariableFeatures(B_TTR_seurat, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(B_PTR_seurat, B_TTR_seurat), dims = 1:20)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)

immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

DimPlot(immune.combined, reduction = "umap", split.by = "stim")

FeaturePlot(immune.combined, features = c("221037"), split.by = "stim", max.cutoff = 3, 
            cols = c("grey", "red"),label = TRUE)
VlnPlot(immune.combined, features = c("221037"), split.by = "stim", 
        pt.size = 0, combine = FALSE)

DefaultAssay(immune.combined) <- "RNA"
clueter0.markers <- FindConservedMarkers(immune.combined, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
head(clueter0.markers)
clueter1.markers <- FindConservedMarkers(immune.combined, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
clueter2.markers <- FindConservedMarkers(immune.combined, ident.1 = 2, grouping.var = "stim", verbose = FALSE)
clueter3.markers <- FindConservedMarkers(immune.combined, ident.1 = 3, grouping.var = "stim", verbose = FALSE)
clueter4.markers <- FindConservedMarkers(immune.combined, ident.1 = 4, grouping.var = "stim", verbose = FALSE)

write.table(clueter0.markers,file = "c:/Users/Administrator/Desktop/GSE108989/cluster0.conservemarkers.txt",sep = "\t",col.names = TRUE,row.names = TRUE)
write.table(clueter1.markers,file = "c:/Users/Administrator/Desktop/GSE108989/cluster1.conservemarkers.txt",sep = "\t",col.names = TRUE,row.names = TRUE)
write.table(clueter2.markers,file = "c:/Users/Administrator/Desktop/GSE108989/cluster2.conservemarkers.txt",sep = "\t",col.names = TRUE,row.names = TRUE)
write.table(clueter3.markers,file = "c:/Users/Administrator/Desktop/GSE108989/cluster3.conservemarkers.txt",sep = "\t",col.names = TRUE,row.names = TRUE)
write.table(clueter4.markers,file = "c:/Users/Administrator/Desktop/GSE108989/cluster4.conservemarkers.txt",sep = "\t",col.names = TRUE,row.names = TRUE)

immune.combined$celltype.stim <- paste(Idents(immune.combined), immune.combined$stim, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.stim"
cluster0_PTR_TTR_diff <- FindMarkers(immune.combined, ident.1 = "0_STIM", ident.2 = "0_CTRL", verbose = FALSE)
cluster1_PTR_TTR_diff <- FindMarkers(immune.combined, ident.1 = "1_STIM", ident.2 = "1_CTRL", verbose = FALSE)
cluster2_PTR_TTR_diff <- FindMarkers(immune.combined, ident.1 = "2_STIM", ident.2 = "2_CTRL", verbose = FALSE)
cluster3_PTR_TTR_diff <- FindMarkers(immune.combined, ident.1 = "3_STIM", ident.2 = "3_CTRL", verbose = FALSE)
#cluster4_PTR_TTR_diff <- FindMarkers(immune.combined, ident.1 = "4_STIM", ident.2 = "4_CTRL", verbose = FALSE)
write.table(cluster0_PTR_TTR_diff,file = "c:/Users/Administrator/Desktop/GSE108989/cluster0_PTR_TTR_diff.txt",sep = "\t",col.names = TRUE,row.names = TRUE)
write.table(cluster1_PTR_TTR_diff,file = "c:/Users/Administrator/Desktop/GSE108989/cluster1_PTR_TTR_diff.txt",sep = "\t",col.names = TRUE,row.names = TRUE)
write.table(cluster2_PTR_TTR_diff,file = "c:/Users/Administrator/Desktop/GSE108989/cluster2_PTR_TTR_diff.txt",sep = "\t",col.names = TRUE,row.names = TRUE)
write.table(cluster3_PTR_TTR_diff,file = "c:/Users/Administrator/Desktop/GSE108989/cluster3_PTR_TTR_diff.txt",sep = "\t",col.names = TRUE,row.names = TRUE)

library(ggplot2)
library(cowplot)
#theme_set(theme_cowplot())
#cluster0.cells <- subset(immune.combined, idents = 0)
#Idents(cluster0.cells) <- "stim"
#avg.cluster0.cells <- log1p(AverageExpression(cluster0.cells, verbose = FALSE)$RNA)
#avg.cluster0.cells$gene <- rownames(avg.cluster0.cells)

#cluster1.cell <- subset(immune.combined, idents = 1)
#Idents(cluster1.cell) <- "stim"
#avg.cluster1.cell <- log1p(AverageExpression(cluster1.cell, verbose = FALSE)$RNA)
#avg.cluster1.cell$gene <- rownames(avg.cluster1.cell)

#cluster2.cell <- subset(immune.combined, idents = 2)
#Idents(cluster2.cell) <- "stim"
#avg.cluster2.cell <- log1p(AverageExpression(cluster2.cell, verbose = FALSE)$RNA)
#avg.cluster2.cell$gene <- rownames(avg.cluster2.cell)

#cluster3.cell <- subset(immune.combined, idents = 3)
#Idents(cluster3.cell) <- "stim"
#avg.cluster3.cell <- log1p(AverageExpression(cluster3.cell, verbose = FALSE)$RNA)
#avg.cluster3.cell$gene <- rownames(avg.cluster3.cell)

#cluster4.cell <- subset(immune.combined, idents = 4)
#Idents(cluster4.cell) <- "stim"
#avg.cluster4.cell <- log1p(AverageExpression(cluster4.cell, verbose = FALSE)$RNA)
#avg.cluster4.cell$gene <- rownames(avg.cluster4.cell)

#cluster5.cell <- subset(immune.combined, idents = 5)
#Idents(cluster5.cell) <- "stim"
#avg.cluster5.cell <- log1p(AverageExpression(cluster5.cell, verbose = FALSE)$RNA)
#avg.cluster5.cell$gene <- rownames(avg.cluster5.cell)

#cluster6.cell <- subset(immune.combined, idents = 6)
#Idents(cluster6.cell) <- "stim"
#avg.cluster6.cell <- log1p(AverageExpression(cluster6.cell, verbose = FALSE)$RNA)
#avg.cluster6.cell$gene <- rownames(avg.cluster6.cell)

#cluster7.cell <- subset(immune.combined, idents = 7)
#Idents(cluster7.cell) <- "stim"
#avg.cluster7.cell <- log1p(AverageExpression(cluster7.cell, verbose = FALSE)$RNA)
#avg.cluster7.cell$gene <- rownames(avg.cluster7.cell)

#genes.to.label = "221037"
#p1 <- ggplot(avg.cluster0.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("cluster0")
#p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
#p2 <- ggplot(avg.cluster1.cell, aes(CTRL, STIM)) + geom_point() + ggtitle("cluster1")
#p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
#p3 <- ggplot(avg.cluster2.cell, aes(CTRL, STIM)) + geom_point() + ggtitle("cluster2")
#p3 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)
#p4 <- ggplot(avg.cluster3.cell, aes(CTRL, STIM)) + geom_point() + ggtitle("cluster3")
#p4 <- LabelPoints(plot = p4, points = genes.to.label, repel = TRUE)
#p5 <- ggplot(avg.cluster4.cell, aes(CTRL, STIM)) + geom_point() + ggtitle("cluster4")
#p5 <- LabelPoints(plot = p5, points = genes.to.label, repel = TRUE)
#p6 <- ggplot(avg.cluster5.cell, aes(CTRL, STIM)) + geom_point() + ggtitle("cluster5")
#p6 <- LabelPoints(plot = p6, points = genes.to.label, repel = TRUE)
#p7 <- ggplot(avg.cluster6.cell, aes(CTRL, STIM)) + geom_point() + ggtitle("cluster6")
#p7 <- LabelPoints(plot = p7, points = genes.to.label, repel = TRUE)
#p8 <- ggplot(avg.cluster7.cell, aes(CTRL, STIM)) + geom_point() + ggtitle("cluster7")
#p8 <- LabelPoints(plot = p8, points = genes.to.label, repel = TRUE)

#p1
#p2
#p3
#p4
#p5
#p6
#p7
#p8

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

FeaturePlot(B_PTR_TTR_Seurat,features = "221037")
VlnPlot(B_PTR_TTR_Seurat,features = "221037",pt.size = 0)

PTR_TTR.markers <- FindAllMarkers(B_PTR_TTR_Seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
PTR_TTR.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(PTR_TTR.markers,file = "c:/Users/Administrator/Desktop/GSE108989/PTR_TTR.markers.txt",sep = "\t",col.names = TRUE,row.names = TRUE)

top10 <- PTR_TTR.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(B_PTR_TTR_Seurat, features = top10$gene) + NoLegend()

new.cluster.ids <- c("control", "TTR", "control", "control", "TTR", "control","control", "control")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
DimPlot(immune.combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
FeaturePlot(immune.combined,features = "221037",label = TRUE)
VlnPlot(immune.combined,features = "221037",pt.size = 0,split.by = "stim")
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:20)
FeaturePlot(immune.combined,features = "221037",reduction = "tsne",label = TRUE)

FeaturePlot(B_PTR_TTR_Seurat,features = "221037",cols = c("yellow","purple"),min.cutoff = "q40",max.cutoff = "q90")

#K<-as.data.frame(B_PTR_TTR_Seurat@assays[["RNA"]]@scale.data)
A_PTR_TTR<-A_PTR_TTR[,-1]
A_PTR_TTR<-A_PTR_TTR[,-1]
rownames(A_PTR_TTR)<-rowname_A
#B_PTR_TTR_Seurat@assays[["RNA"]]@scale.data<-as.matrix(A_PTR_TTR)
#FeaturePlot(B_PTR_TTR_Seurat,features = "221037",cols = c("yellow","purple"),min.cutoff = -5,max.cutoff = 5)

A_PTR_TTR_Seurat<-CreateSeuratObject(counts = A_PTR_TTR, project = "Treg", min.cells = 3, min.features = 200)
#A_PTR_TTR_Seurat <- subset(A_PTR_TTR_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000)
#A_PTR_TTR_Seurat <- NormalizeData(A_PTR_TTR_Seurat)
#A_PTR_TTR_Seurat <- FindVariableFeatures(A_PTR_TTR_Seurat, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(A_PTR_TTR_Seurat)
#A_PTR_TTR_Seurat <- ScaleData(A_PTR_TTR_Seurat, features = all.genes)
A_PTR_TTR_Seurat@assays[["RNA"]]@scale.data<-as.matrix(A_PTR_TTR)
#A_PTR_TTR_Seurat <- RunPCA(A_PTR_TTR_Seurat, features = VariableFeatures(object = A_PTR_TTR_Seurat))
#ElbowPlot(A_PTR_TTR_Seurat)
#A_PTR_TTR_Seurat <- FindNeighbors(A_PTR_TTR_Seurat, dims = 1:10)
#A_PTR_TTR_Seurat <- FindClusters(A_PTR_TTR_Seurat, resolution = 0.5)

PTR<-colnames(A_PTR)
TTR<-colnames(A_TTR)

Idents(A_PTR_TTR_Seurat, cells = PTR)<-'PTR'
Idents(A_PTR_TTR_Seurat, cells = TTR)<-'TTR'
A_PTR_TTR_Seurat <- RunPCA(A_PTR_TTR_Seurat, features = rownames(A_PTR_TTR))
A_PTR_TTR_Seurat <- RunUMAP(A_PTR_TTR_Seurat, dims = 1:10)
A_PTR_TTR_Seurat <- RunTSNE(A_PTR_TTR_Seurat, dims = 1:10)

DimPlot(A_PTR_TTR_Seurat, reduction = "umap")
DimPlot(A_PTR_TTR_Seurat, reduction = "tsne")

FeaturePlot(A_PTR_TTR_Seurat,features = "221037",cols = c("yellow","purple"),max.cutoff = "q95")

FeaturePlot(A_PTR_TTR_Seurat,features = "221037",cols = c("#87CEFA","#FFFF00","#FF4500"),max.cutoff = 5,min.cutoff = -1)
FeaturePlot(A_PTR_TTR_Seurat,features = "221037",cols = c("#191970","#FF4500"),max.cutoff = "q95")
VlnPlot(A_PTR_TTR_Seurat,features = "221037",pt.size = 0)
P<-as.data.frame(A_PTR_TTR_Seurat@assays[["RNA"]]@scale.data)

Correlation<-B_TTR
Correlation<-as.data.frame(Correlation)

Correlation_t<-t(Correlation)
Correlation_t<-as.data.frame(Correlation_t)
write.table(Correlation,file = "c:/Users/Administrator/Desktop/GSE108989/Correlation.txt",sep = "\t",row.names = TRUE,col.names = TRUE)
Foxp3_jmjd1c<-read.table("c:/Users/Administrator/Desktop/GSE108989/Foxp3_Jmjd1c.txt",sep = "\t",row.names = 1,header = TRUE)
#R<-rownames(Foxp3_jmjd1c)

library(ggplot2)
library(ggpubr)
ggplot(data = Foxp3_jmjd1c,aes(x=Foxp3,y=Jmjd1c)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = Foxp3_jmjd1c,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
FeaturePlot(pbmc,features = c("Kdm6b","Nr4a1"),cells = R,label = TRUE,pt.size = 3.4)
