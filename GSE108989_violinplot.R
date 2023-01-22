A<-read.table("c:/Users/Administrator/Desktop/GSE108989/GSE108989_CRC.TCell.S10805.norm.centered.txt",header = TRUE,sep = "\t")
B<-read.table("c:/Users/Administrator/Desktop/GSE108989/GSE108989_CRC.TCell.S11138.count.txt",header = TRUE,sep = "\t")
C<-read.table("c:/Users/Administrator/Desktop/GSE108989/GSE108989_CRC.TCell.S11138.TPM.txt",header = TRUE,sep = "\t")

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


A_PTR_TTR_JMJD1C<-filter(A_PTR_TTR,geneSymbol=="JMJD1C")
B_PTR_TTR_JMJD1C<-filter(B_PTR_TTR,symbol=="JMJD1C")
C_PTR_TTR_JMJD1C<-filter(C_PTR_TTR,symbol=="JMJD1C")

A_PTR_TTR_JMJD1C<-A_PTR_TTR_JMJD1C[,-1]
A_PTR_TTR_JMJD1C<-A_PTR_TTR_JMJD1C[,-1]
A_PTR_TTR_JMJD1C_t<-as.data.frame(t(A_PTR_TTR_JMJD1C))



B_PTR_TTR_JMJD1C<-B_PTR_TTR_JMJD1C[,-1]
B_PTR_TTR_JMJD1C<-B_PTR_TTR_JMJD1C[,-1]
B_PTR_TTR_JMJD1C_t<-as.data.frame(t(B_PTR_TTR_JMJD1C))
#B_PTR_TTR_JMJD1C_t<-10*sqrt(B_PTR_TTR_JMJD1C_t+1000)

C_PTR_TTR_JMJD1C<-C_PTR_TTR_JMJD1C[,-1]
C_PTR_TTR_JMJD1C<-C_PTR_TTR_JMJD1C[,-1]
C_PTR_TTR_JMJD1C_t<-as.data.frame(t(C_PTR_TTR_JMJD1C))

A_PTR_TTR_Group<-c(rep("PTR",836),rep("TTR",1265))
B_PTR_TTR_Group<-c(rep("PTR",844),rep("TTR",1300))
C_PTR_TTR_Group<-c(rep("PTR",844),rep("TTR",1300))

A_PTR_TTR_plot<-data.frame(data=A_PTR_TTR_JMJD1C_t[,1],group=A_PTR_TTR_Group)
B_PTR_TTR_plot<-data.frame(data=B_PTR_TTR_JMJD1C_t[,1],group=B_PTR_TTR_Group)
C_PTR_TTR_plot<-data.frame(data=C_PTR_TTR_JMJD1C_t[,1],group=C_PTR_TTR_Group)

color_A<-factor(A_PTR_TTR_Group)
color_B<-factor(B_PTR_TTR_Group)
color_C<-factor(C_PTR_TTR_Group)

library(ggplot2)
median.quartile <- function(x){
  out <- quantile(x, probs = 0.5)
  names(out) <- "y"
  return(out) 
}
P_A<- ggplot(A_PTR_TTR_plot, aes(x=group, y=data,fill=color_A)) + 
  geom_violin(trim=FALSE,color="blue") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  #geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("") + stat_summary(fun=median.quartile,geom = "")
P_A

P_B<- ggplot(B_PTR_TTR_plot, aes(x=group, y=data,fill=color_B)) + 
  geom_violin(trim=FALSE,color="blue") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("") #设置x轴和y轴的标题
P_B

P_C<- ggplot(C_PTR_TTR_plot, aes(x=group, y=data,fill=color_C)) + 
  geom_violin(trim=FALSE,color="blue") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(angle=15,hjust = 1,colour="black",family="Times",size=20), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=16,face="plain"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 20,face="plain"), #设置y轴标题的字体属性
        panel.border = element_blank(),axis.line = element_line(colour = "black",size=1), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="italic", family="Times", colour="black",  #设置图例的子标题的字体属性
                                 size=16),
        legend.title=element_text(face="italic", family="Times", colour="black", #设置图例的总标题的字体属性
                                  size=18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Value")+xlab("") #设置x轴和y轴的标题
P_C

