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

#------------------------------------------------------------
#    Supplementary Figure 2 NO.1
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
          x = "CNA_loss_Nwgd", 
          y = "CNA_loss_WGD",
          color = "red",
          size = 1.3,)+  
  theme_par()+
  geom_abline(intercept=0,slope=1, cex = 0.8,color = "green")+
  geom_text_repel(aes(CNA_loss_Nwgd, 
                      CNA_loss_WGD,
                      label=rownames(CN_Chrom_mean)),
                  max.overlaps = getOption("ggrepel.max.overlaps", 
                                           default = 100))+
  ggtitle("CNV From Deletion")+
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
#    Supplementary Figure 2 NO.2
#------------------------------------------------------------

ggscatter(CN_Chrom_mean, 
          x = "CNA_gain_nWGD", 
          y = "CNA_gain_WGD",
          color = "red",
          size = 1.3)+  
  theme_par()+
  geom_abline(intercept=0,slope=1, cex = 0.8,color = "green")+
  geom_text_repel(aes(CNA_gain_nWGD, 
                      CNA_gain_WGD,
                      label=rownames(CN_Chrom_mean)),
                  max.overlaps = getOption("ggrepel.max.overlaps",default = 100))+
  ggtitle("CNV From Amplification")+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15),
        legend.position = "none")+
  labs(y="WGD+",x="WGD-")
