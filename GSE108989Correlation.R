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

B_PTR_seurat  <- CreateSeuratObject(counts =B_PTR, project = "B_PTR", min.cells = 3, min.features = 200)
VlnPlot(B_PTR_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
B_PTR_seurat <- subset(B_PTR_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 )
B_PTR_seurat <- NormalizeData(B_PTR_seurat, verbose = FALSE)
B_PTR_seurat <- FindVariableFeatures(B_PTR_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(B_PTR_seurat)
B_PTR_seurat <- ScaleData(B_PTR_seurat, features = all.genes)
B_PTR_seurat <- RunPCA(B_PTR_seurat, features = VariableFeatures(object = B_PTR_seurat))
ElbowPlot(B_PTR_seurat)
B_PTR_seurat <- FindNeighbors(B_PTR_seurat, dims = 1:10)
B_PTR_seurat <- FindClusters(B_PTR_seurat, resolution = 0.5)
B_PTR_seurat <- RunUMAP(B_PTR_seurat, dims = 1:10)
DimPlot(B_PTR_seurat, reduction = "umap")

B_PTR_nomrlize_scaledata<-B_PTR_seurat@assays[["RNA"]]@scale.data
B_PTR_nomrlize_scaledata<-as.data.frame(B_PTR_nomrlize_scaledata)

B_PTR_nomrlize_scaledata_t<-t(B_PTR_nomrlize_scaledata)
B_PTR_nomrlize_scaledata_t<-as.data.frame(B_PTR_nomrlize_scaledata_t)

IL2_stat5<-read.table("g:/GSE98638/IL2_STAT5_entrezID.txt",sep = "\t",header = TRUE)
IL6_stat3<-read.table("g:/GSE98638/IL6_STAT3_entrezID.txt",sep = "\t",header = TRUE)

IL2_stat5<-IL2_stat5[2:200,]
IL6_stat3<-IL6_stat3[2:88,]

B_PTR_seurat <- CellCycleScoring(B_PTR_seurat, s.features = IL2_stat5, g2m.features = IL6_stat3, set.ident = FALSE)
IL2_stat5_list<-list(IL2_stat5)
IL6_stat3_list<-list(IL6_stat3)
B_PTR_seurat <- AddModuleScore(B_PTR_seurat,features = IL2_stat5_list,name = "IL2_stat5")
B_PTR_seurat <- AddModuleScore(B_PTR_seurat,features = IL6_stat3_list,name = "IL6_stat3")
JMJD1C<-select(B_PTR_nomrlize_scaledata_t,starts_with("221037"))
NRP1<-select(B_PTR_nomrlize_scaledata_t,starts_with("8829"))
Pdcd1<-select(B_PTR_nomrlize_scaledata_t,starts_with("5133") & ends_with("5133"))
Ifng<-select(B_PTR_nomrlize_scaledata_t,starts_with("3458") & ends_with("3458"))
IL2_stat5_score<-as.data.frame(B_PTR_seurat$IL2_stat51)
IL6_stat3_score<-as.data.frame(B_PTR_seurat$IL6_stat31)
data<-cbind(IL2_stat5_score,IL6_stat3_score,JMJD1C,Pdcd1,Ifng,NRP1)
colnames(data)<-c("IL2_stat5","IL6_stat3","JMJD1C","Pdcd1","Ifng","NRP1")
#library(corrplot)
#corrplot(matrix, method = "color")  

library(ggplot2)
library(ggpubr)
ggplot(data = data,aes(x=IL2_stat5,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new<-filter(data,IL2_stat5>0 & JMJD1C>0)
ggplot(data = data_new,aes(x=IL2_stat5,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=IL6_stat3,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new2<-filter(data,IL6_stat3>0 & JMJD1C>0)
ggplot(data = data_new2,aes(x=IL6_stat3,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new2,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=Pdcd1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new3<-filter(data,Pdcd1>0 & JMJD1C>0)
ggplot(data = data_new3,aes(x=Pdcd1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new3,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=Ifng,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new4<-filter(data,Ifng>0 & JMJD1C>0)
ggplot(data = data_new4,aes(x=Ifng,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new4,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=NRP1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new5<-filter(data,NRP1>0 & JMJD1C>0)
ggplot(data = data_new5,aes(x=NRP1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new5,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 


B_TTR_seurat  <- CreateSeuratObject(counts =B_TTR, project = "B_TTR", min.cells = 3, min.features = 200)
VlnPlot(B_TTR_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
B_TTR_seurat <- subset(B_TTR_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 )
B_TTR_seurat <- NormalizeData(B_TTR_seurat, verbose = FALSE)
B_TTR_seurat <- FindVariableFeatures(B_TTR_seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(B_TTR_seurat)
B_TTR_seurat <- ScaleData(B_TTR_seurat, features = all.genes)
B_TTR_seurat <- RunPCA(B_TTR_seurat, features = VariableFeatures(object = B_TTR_seurat))
ElbowPlot(B_TTR_seurat)
B_TTR_seurat <- FindNeighbors(B_TTR_seurat, dims = 1:10)
B_TTR_seurat <- FindClusters(B_TTR_seurat, resolution = 0.5)
B_TTR_seurat <- RunUMAP(B_TTR_seurat, dims = 1:10)
DimPlot(B_TTR_seurat, reduction = "umap")

B_TTR_nomrlize_scaledata<-B_TTR_seurat@assays[["RNA"]]@scale.data
B_TTR_nomrlize_scaledata<-as.data.frame(B_TTR_nomrlize_scaledata)

B_TTR_nomrlize_scaledata_t<-t(B_TTR_nomrlize_scaledata)
B_TTR_nomrlize_scaledata_t<-as.data.frame(B_TTR_nomrlize_scaledata_t)

IL2_stat5<-read.table("g:/GSE98638/IL2_STAT5_entrezID.txt",sep = "\t",header = TRUE)
IL6_stat3<-read.table("g:/GSE98638/IL6_STAT3_entrezID.txt",sep = "\t",header = TRUE)

IL2_stat5<-IL2_stat5[2:200,]
IL6_stat3<-IL6_stat3[2:88,]

B_TTR_seurat <- CellCycleScoring(B_TTR_seurat, s.features = IL2_stat5, g2m.features = IL6_stat3, set.ident = FALSE)
IL2_stat5_list<-list(IL2_stat5)
IL6_stat3_list<-list(IL6_stat3)
B_TTR_seurat <- AddModuleScore(B_TTR_seurat,features = IL2_stat5_list,name = "IL2_stat5")
B_TTR_seurat <- AddModuleScore(B_TTR_seurat,features = IL6_stat3_list,name = "IL6_stat3")
JMJD1C<-select(B_TTR_nomrlize_scaledata_t,starts_with("221037"))
NRP1<-select(B_TTR_nomrlize_scaledata_t,starts_with("8829"))
Pdcd1<-select(B_TTR_nomrlize_scaledata_t,starts_with("5133") & ends_with("5133"))
Ifng<-select(B_TTR_nomrlize_scaledata_t,starts_with("3458") & ends_with("3458"))
IL2_stat5_score<-as.data.frame(B_TTR_seurat$IL2_stat51)
IL6_stat3_score<-as.data.frame(B_TTR_seurat$IL6_stat31)
data<-cbind(IL2_stat5_score,IL6_stat3_score,JMJD1C,Pdcd1,Ifng,NRP1)
colnames(data)<-c("IL2_stat5","IL6_stat3","JMJD1C","Pdcd1","Ifng","NRP1")
#library(corrplot)
#corrplot(matrix, method = "color")  

library(ggplot2)
library(ggpubr)
ggplot(data = data,aes(x=IL2_stat5,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new<-filter(data,IL2_stat5>0 & JMJD1C>0)
ggplot(data = data_new,aes(x=IL2_stat5,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=IL6_stat3,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new2<-filter(data,IL6_stat3>0 & JMJD1C>0)
ggplot(data = data_new2,aes(x=IL6_stat3,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new2,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=Pdcd1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new3<-filter(data,Pdcd1>0 & JMJD1C>0)
ggplot(data = data_new3,aes(x=Pdcd1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new3,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=Ifng,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new4<-filter(data,Ifng>0 & JMJD1C>0)
ggplot(data = data_new4,aes(x=Ifng,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new4,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
ggplot(data = data,aes(x=NRP1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 
data_new5<-filter(data,NRP1>0 & JMJD1C>0)
ggplot(data = data_new5,aes(x=NRP1,y=JMJD1C)) + geom_point(color="blue",size=3.4) + stat_smooth(method = "lm",se=FALSE,size=1.5,color="red")+stat_cor(data = data_new5,method = "pearson") +  theme_set(theme_bw()) +  theme_set(theme_bw()) 

