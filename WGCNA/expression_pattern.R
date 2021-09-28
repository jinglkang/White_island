library(ggplot2)
library(gridExtra)
library(ggsci)
library("scales")
library(SummarizedExperiment)
library(DEGreport)
show_col(pal_aaas("default", alpha = 1)(10))
setwd("~/Documents/2021/White_island/WGCNA/Common/")
dat1<-read.table("geneInfo_pH_sig_pos_plot.txt",header=TRUE)
names(dat1)
ggplot(dat1,aes(pH,Nb)) +
  geom_boxplot(aes(group=pH,colour=Site), coef = 6, outlier.shape=NA) +
  geom_jitter(aes(colour=Site),size=0.5,width = 0.002) + 
  #  stat_summary(fun=median, geom="line", aes(group=1))+
  #  geom_line(aes(colour=Site)) +
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(colour="black",family="Times",size=13), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=13,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="plain"), #设置y轴标题的字体属性
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Times", colour="black", size=13),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black", size=13), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Expression Ratio")+xlab("")+ #设置x轴和y轴的标题
  theme(legend.position = "top", legend.title = element_blank()) + 
  scale_color_aaas() + geom_smooth(method=lm, fullrange=TRUE, colour="#8081807F", se=F) +
  annotate(x=8, y=0.15, label=paste("R = ", round(cor(dat1$pH, dat1$Nb),2),
                                       ", p value = ", signif(cor.test(dat1$pH, dat1$Nb)$p.value, 3)), 
           geom="text", size=5, colour="#1B1919FF", family="Times")

dat1<-read.table("geneInfo_pH_sig_neg_plot.txt",header=TRUE)
names(dat1)
ggplot(dat1,aes(pH,Nb)) +
  #  geom_point(aes(colour=Site),size=1) + 
  geom_boxplot(aes(group=pH,colour=Site), coef = 6, outlier.shape=NA) +
  geom_jitter(aes(colour=Site),size=0.5,width = 0.002) + 
  #  stat_summary(fun=median, geom="line", aes(group=1))+
  #  geom_line(aes(colour=Site)) +
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(colour="black",family="Times",size=13), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=13,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="plain"), #设置y轴标题的字体属性
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Times", colour="black", size=10),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Times", colour="black", size=13), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("Expression Ratio")+xlab("")+ #设置x轴和y轴的标题
  theme(legend.position = "top", legend.title = element_blank()) + 
  scale_color_aaas() + geom_smooth(method=lm, fullrange=TRUE, colour="#8081807F") +
  annotate(x=8, y=0.08, label=paste("R = ", round(cor(dat1$pH, dat1$Nb),2),
                                      ", p value = ", signif(cor.test(dat1$pH, dat1$Nb)$p.value, 3)), 
           geom="text", size=5, colour="#1B1919FF", family="Times")
