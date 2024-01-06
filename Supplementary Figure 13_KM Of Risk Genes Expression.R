######################################## KM
library(survminer)
library(survival)
library(rms)
library(ggrisk)
#################################   training
load("training.RData")
genes<-c("ANLN","BIRC3","C9orf68","CAMK4","CD200R1",
         "CD226","CD79A","CLEC10A","DTHD1","EOMES",
         "FCER2","FGL2","KLRB1","LAT","LGALS2",
         "NCR3","PPP2R2B","PSTPIP1","PTPRC","SLA2","SPN","STAP1"  )


splots <- list()
for (i in 1:22) {
  col<-c(genes[i],"OS","OS.time")
  dat3<-dat2[,col]
  dat3$group[dat3[,1] < median(dat3[,1])] <- "Low"
  dat3$group[dat3[,1] >= median(dat3[,1])] <- "High"
  
  fit<-survfit(Surv(OS.time,OS)~group,data=dat3)
  p<-ggsurvplot(fit,
                conf.int = F,# Show confidence intervals          
                title = colnames(dat3)[1], # Setting the title         
                linetype = "strata", # Automatic setting of curve types according to WGD grouping
                palette=c("#EE2C2C","#27408B"), 
                xlab = "Days",# Specify x-axis labels          
                legend = c(0.8,0.9), # Specify legend location    
                legend.title = "", # Setting the legend title  
                ggtheme = theme(axis.line =  element_line(colour = "black"),
                                panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank(),
                                panel.background = element_rect(fill = NA,
                                                                colour = NA),
                                legend.background = element_rect(fill = NA),
                                axis.text = element_text(size=18),
                                axis.title = element_text(size=18),
                                title = element_text(size=21),
                                legend.text = element_text(size=21),
                                legend.key = element_rect(color = NA, # filigree color
                                                          fill = NA), # filler color
                                legend.key.size = unit(0.4, "inches"),
                                plot.margin = margin(1,1,1,1, unit = "cm")
                ), # Setting up the ggplot2 theme
                legend.labs = c("High exp", "Low exp"),
                pval = TRUE,
                pval.size = 9)
  assign(paste("p",i,sep=""),p)
}


splots[[1]] <- p1
splots[[6]] <- p2
splots[[11]] <- p3
splots[[16]] <- p4

splots[[2]] <- p5
splots[[7]] <- p6
splots[[12]] <- p7
splots[[17]] <- p8

splots[[3]] <- p9
splots[[8]] <- p10
splots[[13]] <- p11
splots[[18]] <- p12

splots[[4]] <- p13
splots[[9]] <- p14
splots[[14]] <- p15
splots[[19]] <- p16

splots[[5]] <- p17
splots[[10]] <- p18
splots[[15]] <- p19
splots[[20]] <- p20

splots[[21]] <- p21
splots[[22]] <- p22

#------------------------------------------------------------
#    Supplementary Figure 13 
#------------------------------------------------------------

arrange_ggsurvplots(splots, ncol = 5, nrow = 5) 

#-----------------------------------------------------------
#-----------------------------------------------------------








