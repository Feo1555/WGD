rm(list=ls())

getwd()
setwd("")
library(ggrepel)
library(ggpubr)
library(openxlsx)
library(ggplot2)
library(tidybayes)
library(ggthemes)
library(ggsci)
library(tidyverse)
library(plyr)

#------------------------------------------------------------
#    Figure 1A NO.1
#------------------------------------------------------------

arm<-read.xlsx("Supplementary File-data1.xlsx",sheet = 1)
rownames(arm)<-arm[,1]
colnames(arm)<-arm[1,]
arm<-arm[-1,-1]
Data<-as.data.frame(t(arm))
colnames(Data)[1]<-"WGD"
Data$burden<-as.numeric(Data$burden)
Data$WGD[Data$WGD %in% "0"] <- "WGD-"
Data$WGD[Data$WGD %in% "1"] <- "WGD+"
Data$WGD<-factor(Data$WGD, levels=c('WGD+', 'WGD-'))
Data_summary <- summarySE(Data, measurevar="burden", groupvars="WGD")

ggplot(Data, aes(x=WGD, y=burden,fill=WGD)) + 
  stat_eye(aes(alpha = stat(f)), 
           show_interval = FALSE, 
           trim=FALSE, 
           position=position_dodge(0.9)) + 
  geom_point(data = Data_summary,
             aes(x=WGD, 
                 y=burden),
             pch=19,
             position=position_dodge(0.9),
             size=1.5)+ 
  geom_errorbar(data = Data_summary,
                aes(ymin = burden-sd,
                    ymax = burden+sd), 
                width=0.1,
                position=position_dodge(0.9),
                color="black",
                alpha = 0.7,
                size=0.5) +
  labs(x="",y="Aneuploidy Burden on Arm")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 12)+ 
  scale_fill_manual(values = c("#EE3B3B","#3A5FCD"))+ 
  theme_par()+
  #facet_wrap(.~cpd_name,scales = "free_y",nrow=1)+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 35),
        legend.position = 'none',
        panel.border = element_blank(),
        #legend.title=element_blank(),
        #legend.text = element_text(size=10),
        #axis.text = element_text(size=15),
        #axis.title = element_text(size=15),
        #legend.key.size = unit(0.3, "inches")
  )

#------------------------------------------------------------
#    Figure 1A NO.2
#------------------------------------------------------------

a<-read.xlsx("Supplementary File-data1.xlsx",sheet = 2)
rownames(a)<-a[,1]
colnames(a)<-a[1,]
a<-a[-1,-1]
chrom<-rownames(a)
a[,1]<-as.numeric(a[,1])
a[,2]<-as.numeric(a[,2])

ggscatter(a, 
          x = "nWGD", 
          y = "WGD",
          color = "red",
          size = 1.3,)+  
  theme_par()+
  geom_abline(intercept=0,slope=1, cex = 0.8,color = "gray")+
  geom_text_repel(aes(nWGD, 
                      WGD,
                      label=rownames(a)),
                  max.overlaps = getOption("ggrepel.max.overlaps", 
                                           default = 100))+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15),
        legend.position = "none")+
  labs(y="WGD+",x="WGD-")

#------------------------------------------------------------
#    Figure 1A NO.3
#------------------------------------------------------------

chrom<-read.xlsx("Supplementary File-data1.xlsx",sheet = 3)
rownames(chrom)<-chrom[,1]
colnames(chrom)<-chrom[1,]
chrom<-chrom[-1,-1]
Data<-as.data.frame(t(chrom))
colnames(Data)[1]<-"WGD"
Data$burden<-as.numeric(Data$burden)
Data$WGD[Data$WGD %in% "0"] <- "WGD-"
Data$WGD[Data$WGD %in% "1"] <- "WGD+"
Data$WGD<-factor(Data$WGD, levels=c('WGD+', 'WGD-'))
Data_summary <- summarySE(Data, measurevar="burden", groupvars="WGD")

ggplot(Data, aes(x=WGD, y=burden,fill=WGD)) + 
  stat_eye(aes(alpha = stat(f)), 
           show_interval = FALSE, 
           trim=FALSE, 
           position=position_dodge(0.9)) + 
  geom_point(data = Data_summary,
             aes(x=WGD, 
                 y=burden),
             pch=19,
             position=position_dodge(0.9),
             size=1.5)+ 
  geom_errorbar(data = Data_summary,
                aes(ymin = burden-sd,
                    ymax = burden+sd), 
                width=0.1, 
                position=position_dodge(0.9),
                color="black",
                alpha = 0.7,
                size=0.5) +
  labs(x="",y="Aneuploidy Burden on Chrom")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 12)+ 
  scale_fill_manual(values = c("#EE3B3B","#3A5FCD"))+ 
  theme_par()+
  #facet_wrap(.~cpd_name,scales = "free_y",nrow=1)+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 35),
        legend.position = 'none',
        panel.border = element_blank(),
        #legend.title=element_blank(),
        #legend.text = element_text(size=10),
        #axis.text = element_text(size=15),
        #axis.title = element_text(size=15),
        #legend.key.size = unit(0.3, "inches")
  )

#------------------------------------------------------------
#    Figure 1A NO.4
#------------------------------------------------------------

a<-read.xlsx("Supplementary File-data1.xlsx",sheet = 4)
rownames(a)<-a[,1]
colnames(a)<-a[1,]
a<-a[-1,-1]
chrom<-rownames(a)
a[,1]<-as.numeric(a[,1])
a[,2]<-as.numeric(a[,2])

ggscatter(a, 
          x = "nWGD", 
          y = "WGD",
          color = "red",
          size = 1.3,)+  
  theme_par()+
  geom_abline(intercept=0,slope=1, cex = 0.8,color = "gray")+
  geom_text_repel(aes(nWGD, 
                      WGD,
                      label=rownames(a)),
                  max.overlaps = getOption("ggrepel.max.overlaps", 
                                           default = 100))+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15),
        legend.position = "none")+
  labs(y="WGD+",x="WGD-")

#------------------------------------------------------------
#    Figure 1B NO.1
#------------------------------------------------------------

CIN0<-read.xlsx("Supplementary File-data2.xlsx",sheet = 1)
rownames(CIN0)<-CIN0[,1]
colnames(CIN0)<-CIN0[1,]
CIN0<-CIN0[-1,-1]
CIN0[,4]<-rep("0",nrow(CIN0))
colnames(CIN0)[4]<-"WGD"
colnames(CIN0)[3]<-"CIN"

CIN1<-read.xlsx("Supplementary File-data2.xlsx",sheet = 2)
rownames(CIN1)<-CIN1[,1]
colnames(CIN1)<-CIN1[1,]
CIN1<-CIN1[-1,-1]
CIN1[,4]<-rep("1",nrow(CIN1))
colnames(CIN1)[4]<-"WGD"
colnames(CIN1)[3]<-"CIN"
CIN_WGD<-rbind(CIN1[,c(3,4)],CIN0[,c(3,4)])
BRCA_CIN<-CIN_WGD
BRCA_CIN[,1]<-as.numeric(BRCA_CIN[,1])


BRCA_CIN$WGD[BRCA_CIN$WGD %in% 0]<-"WGD-"
BRCA_CIN$WGD[BRCA_CIN$WGD %in% 1]<-"WGD+"
BRCA_CIN$WGD<-factor(BRCA_CIN$WGD, levels=c('WGD+', 'WGD-'))



Data_summary <- summarySE(BRCA_CIN, measurevar="CIN", groupvars="WGD")


ggplot(BRCA_CIN, aes(x=WGD, y=CIN,fill=WGD)) + 
  stat_eye(aes(alpha = stat(f)), 
           show_interval = FALSE, 
           trim=FALSE, 
           position=position_dodge(0.9)) + 
  geom_point(data = Data_summary,
             aes(x=WGD, 
                 y=CIN),
             pch=19,
             position=position_dodge(0.9),
             size=1.5)+ 
  geom_errorbar(data = Data_summary,
                aes(ymin = CIN-sd,
                    ymax = CIN+sd), 
                width=0.1, 
                position=position_dodge(0.9),
                color="black",
                alpha = 0.7,
                size=0.5) +
  labs(x="",y="Chromosome Instability")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 12,
                     label.y = 24)+ 
  scale_fill_manual(values = c("#EE3B3B","#3A5FCD"))+ 
  theme_par()+
  #facet_wrap(.~cpd_name,scales = "free_y",nrow=1)+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 35),
        legend.position = 'none',
        panel.border = element_blank(),
        #legend.title=element_blank(),
        #legend.text = element_text(size=10),
        #axis.text = element_text(size=15),
        #axis.title = element_text(size=15),
        #legend.key.size = unit(0.3, "inches")
  )

#------------------------------------------------------------
#    Figure 1B NO.2
#------------------------------------------------------------

a<-read.csv("CIN_Chrom_mean.csv",row.names = 1)
data<-a[,c(3,6)]
data$WGD1<-"WGD-"
data$WGD2<-"WGD+"
data1<-data[,c(1,3)]
data1$chrom<-rownames(data1)
rownames(data1)<-1:22
data2<-data[,c(2,4)]
data2$chrom<-rownames(data2)
rownames(data2)<-1:22
colnames(data1)[1:2]<-c("BRCA_CIN","WGD")
colnames(data2)[1:2]<-c("BRCA_CIN","WGD")
Data<-rbind(data2,data1)
Data$WGD<-factor(Data$WGD,levels = c("WGD+","WGD-"))



Data_summary <- summarySE(Data, measurevar="BRCA_CIN", groupvars="WGD")
ggplot(Data, aes(x=WGD, y=BRCA_CIN,fill=WGD)) + 
  stat_eye(aes(alpha = stat(f)), 
           show_interval = FALSE, 
           trim=FALSE, 
           side = "left",
           position=position_dodge(0.9)) + 
  geom_point(data = Data_summary,
             aes(x=WGD, 
                 y=BRCA_CIN),
             pch=19,
             position=position_dodge(0.9),
             size=1.5)+ 
  geom_errorbar(data = Data_summary,
                aes(ymin = BRCA_CIN-sd,
                    ymax = BRCA_CIN+sd), 
                width=0, 
                position=position_dodge(0.9),
                color="black",
                alpha = 0.7,
                size=0.5) +
  labs(x="",y="CIN on Chrom")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 12,
                     label.y = 0.98)+ 
  geom_line(aes(group = chrom,color = chrom), lwd = 0.6)+
  scale_fill_manual(values = c("#EE3B3B","#3A5FCD"))+ 
  theme_par()+
  #facet_wrap(.~cpd_name,scales = "free_y",nrow=1)+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 35),
        legend.position = 'none',
        panel.border = element_blank(),
        #legend.title=element_blank(),
        #legend.text = element_text(size=10),
        #axis.text = element_text(size=15),
        #axis.title = element_text(size=15),
        #legend.key.size = unit(0.3, "inches")
  )

#------------------------------------------------------------
#    Figure 1C NO.1
#------------------------------------------------------------

library(openxlsx)
CNV0<-read.xlsx("CN_WGD0.xlsx",rowNames = T,sheet = "BRCA")
CNV0<-CNV0[,c(4,5,6)]
CNV0$WGD<-rep("0",nrow(CNV0))
CNV0$seg_total<-CNV0$`num(loss_seg)`+CNV0$`num(gain_seg)`
colnames(CNV0)[1:3]<-c("CNV_loss","CNV_gain","CNV_total")

CNV1<-read.xlsx("CN_WGD1.xlsx",rowNames = T,sheet = "BRCA")
CNV1<-CNV1[,c(4,5,6)]
CNV1$WGD<-rep("1",nrow(CNV1))
CNV1$seg_total<-CNV1$`num(loss_seg)`+CNV1$`num(gain_seg)`
colnames(CNV1)[1:3]<-c("CNV_loss","CNV_gain","CNV_total")

CNV_WGD<-rbind(CNV1,CNV0)
CNV_WGD[,1:3]<-log10(CNV_WGD[,1:3])

CNV_WGD$CNV_loss[CNV_WGD$CNV_loss %in% -Inf] <- 0
CNV_WGD$CNV_gain[CNV_WGD$CNV_gain %in% -Inf] <- 0
CNV_WGD$CNV_total[CNV_WGD$CNV_total %in% -Inf] <- 0

CNV_WGD$WGD[CNV_WGD$WGD == 0]<-"WGD-"
CNV_WGD$WGD[CNV_WGD$WGD == 1]<-"WGD+"
CNV_WGD$WGD<-factor(CNV_WGD$WGD, levels=c('WGD+', 'WGD-'))

Data_summary <- summarySE(CNV_WGD, measurevar="CNV_total", groupvars="WGD")
ggplot(CNV_WGD, aes(x=WGD, y=CNV_total,fill=WGD)) + 
  stat_eye(aes(alpha = stat(f)), 
           show_interval = FALSE, 
           trim=FALSE, 
           position=position_dodge(0.9)) + 
  geom_point(data = Data_summary,
             aes(x=WGD, 
                 y=CNV_total),
             pch=19,
             position=position_dodge(0.9),
             size=1.5)+ 
  geom_errorbar(data = Data_summary,
                aes(ymin = CNV_total-sd,
                    ymax = CNV_total+sd), 
                width=0.1, 
                position=position_dodge(0.9),
                color="black",
                alpha = 0.7,
                size=0.5) +
  labs(x="",y="Log10(CNV)")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 12)+ 
  scale_fill_manual(values = c("#EE3B3B","#3A5FCD"))+ 
  theme_par()+
  #facet_wrap(.~cpd_name,scales = "free_y",nrow=1)+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 35),
        legend.position = 'none',
        panel.border = element_blank(),
        #legend.title=element_blank(),
        #legend.text = element_text(size=10),
        #axis.text = element_text(size=15),
        #axis.title = element_text(size=15),
        #legend.key.size = unit(0.3, "inches")
  )

#------------------------------------------------------------
#    Figure 1C NO.2
#------------------------------------------------------------

CN_Chrom_mean<-read.xlsx("Supplementary File-data3.xlsx",sheet = 3)
rownames(CN_Chrom_mean)<-CN_Chrom_mean[,1]
colnames(CN_Chrom_mean)<-CN_Chrom_mean[1,]
CN_Chrom_mean<-CN_Chrom_mean[-1,-1]
for (i in seq_along(CN_Chrom_mean)) {
  CN_Chrom_mean[[i]] <- as.numeric(CN_Chrom_mean[[i]])
}
rownames(CN_Chrom_mean)<-gsub("om", "", rownames(CN_Chrom_mean))

ggscatter(CN_Chrom_mean, 
          x = "CNA_both_nWGD", 
          y = "CNA_both_WGD",
          color = "red",
          size = 1.3,)+  
  theme_par()+
  geom_abline(intercept=0,slope=1, cex = 0.8,color = "gray")+
  geom_text_repel(aes(CNA_both_nWGD, 
                      CNA_both_WGD,
                      label=rownames(CN_Chrom_mean)),
                  max.overlaps = getOption("ggrepel.max.overlaps", 
                                           default = 100))+
  ggtitle("")+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15),
        legend.position = "none")+
  labs(y="WGD+",x="WGD-")

#------------------------------------------------------------
#    Figure 1C NO.3
#------------------------------------------------------------

Data_summary <- summarySE(CNV_WGD, measurevar="CNV_gain", groupvars="WGD")
ggplot(CNV_WGD, aes(x=WGD, y=CNV_gain,fill=WGD)) + 
  stat_eye(aes(alpha = stat(f)), 
           show_interval = FALSE, 
           trim=FALSE, 
           position=position_dodge(0.9)) + 
  geom_point(data = Data_summary,
             aes(x=WGD, 
                 y=CNV_gain),
             pch=19,
             position=position_dodge(0.9),
             size=1.5)+ 
  geom_errorbar(data = Data_summary,
                aes(ymin = CNV_gain-sd,
                    ymax = CNV_gain+sd), 
                width=0.1, 
                position=position_dodge(0.9),
                color="black",
                alpha = 0.7,
                size=0.5) +
  labs(x="",y="Log10(CNV from Gain)")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 12,
                     label.y = 2.55)+ 
  scale_fill_manual(values = c("#EE3B3B","#3A5FCD"))+ 
  theme_par()+
  #facet_wrap(.~cpd_name,scales = "free_y",nrow=1)+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 35),
        legend.position = 'none',
        panel.border = element_blank(),
        #legend.title=element_blank(),
        #legend.text = element_text(size=10),
        #axis.text = element_text(size=15),
        #axis.title = element_text(size=15),
        #legend.key.size = unit(0.3, "inches")
  )

#------------------------------------------------------------
#    Figure 1C NO.4
#------------------------------------------------------------

Data_summary <- summarySE(CNV_WGD, measurevar="CNV_loss", groupvars="WGD")
ggplot(CNV_WGD, aes(x=WGD, y=CNV_loss,fill=WGD)) + 
  stat_eye(aes(alpha = stat(f)), 
           show_interval = FALSE, 
           trim=FALSE, 
           position=position_dodge(0.9)) + 
  geom_point(data = Data_summary,
             aes(x=WGD, 
                 y=CNV_loss),
             pch=19,
             position=position_dodge(0.9),
             size=1.5)+ 
  geom_errorbar(data = Data_summary,
                aes(ymin = CNV_loss-sd,
                    ymax = CNV_loss+sd), 
                width=0.1, 
                position=position_dodge(0.9),
                color="black",
                alpha = 0.7,
                size=0.5) +
  labs(x="",y="Log10(CNV from Loss)")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 12)+ 
  scale_fill_manual(values = c("#EE3B3B","#3A5FCD"))+ 
  theme_par()+
  #facet_wrap(.~cpd_name,scales = "free_y",nrow=1)+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 35),
        legend.position = 'none',
        panel.border = element_blank(),
        #legend.title=element_blank(),
        #legend.text = element_text(size=10),
        #axis.text = element_text(size=15),
        #axis.title = element_text(size=15),
        #legend.key.size = unit(0.3, "inches")
  )

#------------------------------------------------------------
#    Figure 1D NO.1
#------------------------------------------------------------



mutation1<-read.csv("adj_mutation-load_updated.csv",row.names = 1)
matr<-subset(mutation1,Cohort %in% "BRCA")
matr$WGD[matr$WGD %in% "nWGD"] <- "WGD-"
matr$WGD[matr$WGD %in% "WGD"] <- "WGD+"
matr$WGD<-factor(matr$WGD, levels=c('WGD+', 'WGD-'))
matr$adj_log10_Non.silent.per.Mb<-matr$log10_Non.silent.per.Mb / matr$ploidy
matr$log10_Non.silent.per.Mb[matr$log10_Non.silent.per.Mb %in% -Inf] <- 0
Data_summary <- summarySE(matr, measurevar="log10_Non.silent.per.Mb", groupvars="WGD")

ggplot(matr, aes(x=WGD, y=log10_Non.silent.per.Mb,fill=WGD)) + 
  stat_eye(aes(alpha = stat(f)), 
           show_interval = FALSE, 
           trim=FALSE, 
           position=position_dodge(0.9)) + 
  geom_point(data = Data_summary,
             aes(x=WGD, 
                 y=log10_Non.silent.per.Mb),
             pch=19,
             position=position_dodge(0.9),
             size=1.5)+ 
  geom_errorbar(data = Data_summary,
                aes(ymin = log10_Non.silent.per.Mb-sd,
                    ymax = log10_Non.silent.per.Mb+sd), 
                width=0.1, 
                position=position_dodge(0.9),
                color="black",
                alpha = 0.7,
                size=0.5) +
  labs(x="",y="Mutation Log10(Non.silent.per.Mb)")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 12)+ 
  scale_fill_manual(values = c("#EE3B3B","#3A5FCD"))+ 
  theme_par()+
  #facet_wrap(.~cpd_name,scales = "free_y",nrow=1)+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 35),
        legend.position = 'none',
        panel.border = element_blank(),
        #legend.title=element_blank(),
        #legend.text = element_text(size=10),
        #axis.text = element_text(size=15),
        #axis.title = element_text(size=15),
        #legend.key.size = unit(0.3, "inches")
  )

#------------------------------------------------------------
#    Figure 1D NO.2
#------------------------------------------------------------

matr$adj_log10_Non.silent.per.Mb[matr$adj_log10_Non.silent.per.Mb %in% -Inf] <- 0
Data_summary <- summarySE(matr, measurevar="adj_log10_Non.silent.per.Mb", groupvars="WGD")

ggplot(matr, aes(x=WGD, y=adj_log10_Non.silent.per.Mb,fill=WGD)) + 
  stat_eye(aes(alpha = stat(f)), 
           show_interval = FALSE, 
           trim=FALSE, 
           position=position_dodge(0.9)) + 
  geom_point(data = Data_summary,
             aes(x=WGD, 
                 y=adj_log10_Non.silent.per.Mb),
             pch=19,
             position=position_dodge(0.9),
             size=1.5)+ 
  geom_errorbar(data = Data_summary,
                aes(ymin = adj_log10_Non.silent.per.Mb-sd,
                    ymax = adj_log10_Non.silent.per.Mb+sd), 
                width=0.1, 
                position=position_dodge(0.9),
                color="black",
                alpha = 0.7,
                size=0.5) +
  labs(x="",y="Mutation adj(Log10(Non.silent.per.Mb))")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 12)+ 
  scale_fill_manual(values = c("#EE3B3B","#3A5FCD"))+
  theme_par()+
  #facet_wrap(.~cpd_name,scales = "free_y",nrow=1)+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 35),
        legend.position = 'none',
        panel.border = element_blank(),
        #legend.title=element_blank(),
        #legend.text = element_text(size=10),
        #axis.text = element_text(size=15),
        #axis.title = element_text(size=15),
        #legend.key.size = unit(0.3, "inches")
  )
