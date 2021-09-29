library(ggplot2)
library(gridExtra)
library(ggsci)
library("scales")
library(SummarizedExperiment)
library(DEGreport)
show_col(pal_aaas("default", alpha = 1)(10))

################################################### 
# Blenny
setwd("~/Documents/2021/White_island/WGCNA/Blenny")
################################################### 
# pH
#### pos
dat1<-read.table("geneInfo_pH_sig_pos_plot.txt",header=TRUE)
names(dat1)
Blenny_pH_pos <- ggplot(dat1,aes(pH,Nb)) +
  #  geom_boxplot(aes(group=pH,colour=Site), coef = 6, outlier.shape=NA, width = 0.02) +
  geom_jitter(aes(colour=Site),size=0.5,width = 0.009, alpha=0.5) + 
  #  stat_summary(fun=median, geom="line", aes(group=1))+
  #  geom_line(aes(colour=Site)) +
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(colour="black",family="Arial",size=13), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Arial大小为20
        axis.text.y=element_text(family="Arial",size=13,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 15,face="plain"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Arial",size = 15,face="plain"),
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black", size=13),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Arial", colour="black", size=13), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("")+xlab("")+ #设置x轴和y轴的标题
  theme(legend.position = "top",legend.title = element_blank()) +
  scale_color_aaas() + geom_smooth(method=lm, fullrange=TRUE, colour="#8081807F", se=F) +
  annotate(x=7.65, y=0.10, label=paste("Blenny: ",length(unique(dat1$Gene)), " genes\n","R = ", round(cor(dat1$pH, dat1$Nb),2),
                                       " (p value = ", signif(cor.test(dat1$pH, dat1$Nb)$p.value, 3),")", sep = ""), 
           geom="text", size=4, colour="#1B1919FF", family="Arial")

Blenny_pH_pos

#### neg
dat1<-read.table("geneInfo_pH_sig_neg_plot.txt",header=TRUE)
names(dat1)
Blenny_pH_neg <- ggplot(dat1,aes(pH,Nb)) +
  #  geom_boxplot(aes(group=pH,colour=Site), coef = 6, outlier.shape=NA, width = 0.02) +
  geom_jitter(aes(colour=Site),size=0.5,width = 0.009, alpha=0.5) + 
  #  stat_summary(fun=median, geom="line", aes(group=1))+
  #  geom_line(aes(colour=Site)) +
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(colour="black",family="Arial",size=13), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Arial大小为20
        axis.text.y=element_text(family="Arial",size=13,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 15,face="plain"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Arial",size = 15,face="plain"),
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black", size=13),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Arial", colour="black", size=13), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("")+xlab("")+ #设置x轴和y轴的标题
  theme(legend.position = "top",legend.title = element_blank()) +
  scale_color_aaas() + geom_smooth(method=lm, fullrange=TRUE, colour="#8081807F", se=F) +
  annotate(x=7.8, y=0.11, label=paste("Blenny: ",length(unique(dat1$Gene)), " genes\n","R = ", round(cor(dat1$pH, dat1$Nb),2),
                                      " (p value = ", signif(cor.test(dat1$pH, dat1$Nb)$p.value, 3),")", sep = ""), 
           geom="text", size=4, colour="#1B1919FF", family="Arial")

Blenny_pH_neg

# Length
#### pos
dat1<-read.table("geneInfo_Length_sig_pos_plot.txt",header=TRUE)
names(dat1)
Blenny_Length_pos <- ggplot(dat1,aes(Length,Nb)) +
  #geom_boxplot(aes(group=Length,colour=Site), coef = 6, outlier.shape=NA, width = 0.15) +
  geom_jitter(aes(colour=Site),size=0.5,width = 0.07, alpha=0.5) + 
  #  stat_summary(fun=median, geom="line", aes(group=1))+
  #  geom_line(aes(colour=Site)) +
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(colour="black",family="Arial",size=13), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Arial大小为20
        axis.text.y=element_text(family="Arial",size=13,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 15,face="plain"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Arial",size = 15,face="plain"),
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black", size=13),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Arial", colour="black", size=13), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("")+xlab("")+ #设置x轴和y轴的标题
  theme(legend.position = "top",legend.title = element_blank()) +
  scale_color_aaas() + geom_smooth(method=lm, fullrange=TRUE, colour="#8081807F", se=F) +
  annotate(x=5, y=0.15, label=paste("Blenny: ",length(unique(dat1$Gene)), " genes\n","R = ", round(cor(dat1$Length, dat1$Nb),2),
                                    " (p value = ", signif(cor.test(dat1$Length, dat1$Nb)$p.value, 3),")", sep = ""), 
           geom="text", size=4, colour="#1B1919FF", family="Arial")

Blenny_Length_pos

#### neg
dat1<-read.table("geneInfo_Length_sig_neg_plot.txt",header=TRUE)
names(dat1)
Blenny_Length_neg <- ggplot(dat1,aes(Length,Nb)) +
  #  geom_boxplot(aes(group=Length,colour=Site), coef = 6, outlier.shape=NA, width = 0.15) +
  geom_jitter(aes(colour=Site),size=0.5,width = 0.07, alpha=0.5) + 
  #  stat_summary(fun=median, geom="line", aes(group=1))+
  #  geom_line(aes(colour=Site)) +
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(colour="black",family="Arial",size=13), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Arial大小为20
        axis.text.y=element_text(family="Arial",size=13,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 15,face="plain"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Arial",size = 15,face="plain"),
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black", size=13),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Arial", colour="black", size=13), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("")+xlab("")+ #设置x轴和y轴的标题
  theme(legend.position = "top",legend.title = element_blank()) +
  scale_color_aaas() + geom_smooth(method=lm, fullrange=TRUE, colour="#8081807F", se=F) +
  annotate(x=6.7, y=0.08, label=paste("Blenny: ",length(unique(dat1$Gene)), " genes\n","R = ", round(cor(dat1$Length, dat1$Nb),2),
                                      " (p value = ", signif(cor.test(dat1$Length, dat1$Nb)$p.value, 3),")", sep = ""), 
           geom="text", size=4, colour="#1B1919FF", family="Arial")

Blenny_Length_neg

# Salinity
#### pos
dat1<-read.table("geneInfo_Salinity_sig_pos_plot.txt",header=TRUE)
names(dat1)
Blenny_Salinity_pos <- ggplot(dat1,aes(Salinity,Nb)) +
  #  geom_boxplot(aes(group=Salinity,colour=Site), coef = 6, outlier.shape=NA, width = 0.15) +
  geom_jitter(aes(colour=Site),size=0.5, width = 0.04, alpha=0.5) + 
  #  stat_summary(fun=median, geom="line", aes(group=1))+
  #  geom_line(aes(colour=Site)) +
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(colour="black",family="Arial",size=13), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Arial大小为20
        axis.text.y=element_text(family="Arial",size=13,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 15,face="plain"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Arial",size = 15,face="plain"),
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black", size=13),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Arial", colour="black", size=13), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("")+xlab("")+ #设置x轴和y轴的标题
  theme(legend.position = "top",legend.title = element_blank()) +
  scale_color_aaas() + geom_smooth(method=lm, fullrange=TRUE, colour="#8081807F", se=F) +
  annotate(x=33.5, y=0.052, label=paste("Blenny: ",length(unique(dat1$Gene)), " genes\n","R = ", round(cor(dat1$Salinity, dat1$Nb),2),
                                        " (p value = ", signif(cor.test(dat1$Salinity, dat1$Nb)$p.value, 3),")", sep = ""), 
           geom="text", size=4, colour="#1B1919FF", family="Arial")

Blenny_Salinity_pos

#### neg
dat1<-read.table("geneInfo_Salinity_sig_neg_plot.txt",header=TRUE)
names(dat1)
Blenny_Salinity_neg <- ggplot(dat1,aes(Salinity,Nb)) +
  #  geom_boxplot(aes(group=Salinity,colour=Site), coef = 6, outlier.shape=NA, width = 0.15) +
  geom_jitter(aes(colour=Site),size=0.5, width = 0.04, alpha=0.5) + 
  #  stat_summary(fun=median, geom="line", aes(group=1))+
  #  geom_line(aes(colour=Site)) +
  theme_bw()+ #背景变为白色
  theme(axis.text.x=element_text(colour="black",family="Arial",size=13), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Arial大小为20
        axis.text.y=element_text(family="Arial",size=13,face="plain",colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Arial",size = 15,face="plain"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Arial",size = 15,face="plain"),
        # panel.border = element_blank(),axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.text=element_text(face="plain", family="Arial", colour="black", size=13),  #设置图例的子标题的字体属性
        legend.title=element_text(face="plain", family="Arial", colour="black", size=13), #设置图例的总标题的字体属性
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank())+  #不显示网格线
  ylab("")+xlab("")+ #设置x轴和y轴的标题
  theme(legend.position = "top",legend.title = element_blank()) +
  scale_color_aaas() + geom_smooth(method=lm, fullrange=TRUE, colour="#8081807F", se=F) +
  annotate(x=33, y=0.2, label=paste("Blenny: ",length(unique(dat1$Gene)), " genes\n","R = ", round(cor(dat1$Salinity, dat1$Nb),2),
                                    " (p value = ", signif(cor.test(dat1$Salinity, dat1$Nb)$p.value, 3),")", sep = ""), 
           geom="text", size=4, colour="#1B1919FF", family="Arial")

Blenny_Salinity_neg
