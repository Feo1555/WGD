library(survminer)
library(survival)
#############################################    Training
#"clinic_train.csv" includes the clinical characteristics of the training group sample
data<-read.csv("clinic_train.csv",row.names = 1,
               check.names = F)
###age
old<-subset(data,Age %in% "Age >= 55")
young<-subset(data,Age %in% "Age < 55")

###stage
data$Stage1[data$Stage %in% "StageI"|
              data$Stage %in% "StageII"] <-"Stage I+II"
data$Stage1[data$Stage %in% "StageIII"|
              data$Stage %in% "StageIV"] <-"Stage III+IV"
table(data$Stage1)
high<-subset(data,Stage1 %in% "Stage III+IV")
low<-subset(data,Stage1 %in% "Stage I+II")
###WGD
WGD<-subset(data,WGD %in% "WGD+")
nWGD<-subset(data,WGD %in% "WGD-")
###subtype
table(data$SubType)
Basal<-subset(data,SubType %in% "Basal")
Her2<-subset(data,SubType %in% "Her2")
LumA<-subset(data,SubType %in% "LumA")
LumB<-subset(data,SubType %in% "LumB")
Normal<-subset(data,SubType %in% "Normal")

#######################################   Age
#------------------------------------------------------------
#    Supplementary Figure 4 NO.1 (left to right )
#------------------------------------------------------------

fit<-survfit(Surv(OS.time,OS)~Group,data=old)
print(fit)
p1<-ggsurvplot(fit,
               conf.int = F,# Show confidence intervals          
               title = "Age > 55", # Setting the title         
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,# Add risk table
               #surv.median.line = "hv", # Setting the Median Survival Display
               xlab = "Days",# Specify x-axis labels          
               legend = c(0.8,0.9), # Specify legend location    
               legend.title = "", # Setting the legend title 
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6), #Set vjust to a negative value
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_rect(fill = "#FFE4E1",
                                                               colour = "#FFE4E1"),
                               legend.background = element_rect(fill = NA),
                               axis.text = element_text(size=15),
                               axis.title = element_text(size=15),
                               title = element_text(size=18),
                               legend.text = element_text(size=18),
                               legend.key = element_rect(color = NA, # filigree color
                                                         fill = NA), # filler color
                               legend.key.size = unit(0.3, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 7)

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Supplementary Figure 4 NO.2
#------------------------------------------------------------

fit<-survfit(Surv(OS.time,OS)~Group,data=young)
print(fit)
p2<-ggsurvplot(fit,
               conf.int = FALSE,# Show confidence intervals         
               title = "Age <= 55", # Setting the title          
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,# Add risk table
               xlab = "Days",# Specify x-axis labels           
               legend = c(0.8,0.9), # Specify legend location    
               legend.title = "", # Setting the legend title        
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_rect(fill = "#FFE4E1",
                                                               colour = "#FFE4E1"),
                               legend.background = element_rect(fill = NA),
                               axis.text = element_text(size=15),
                               axis.title = element_text(size=15),
                               title = element_text(size=18),
                               legend.text = element_text(size=18),
                               legend.key = element_rect(color = NA, # filigree color
                                                         fill = NA), # filler color
                               legend.key.size = unit(0.3, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 7)

#-----------------------------------------------------------
#-----------------------------------------------------------

#######################################   Stage

#------------------------------------------------------------
#    Supplementary Figure 4 NO.3
#------------------------------------------------------------

fit<-survfit(Surv(OS.time,OS)~Group,data=low)
print(fit)
p3<-ggsurvplot(fit,
               conf.int = FALSE,# Show confidence intervals           
               title = "Stage I+II", #  Setting the title          
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,#  Add risk table  
               xlab = "Days",# Specify x-axis labels         
               legend = c(0.8,0.9), # Specify legend location     
               legend.title = "", # Setting the legend title     
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_rect(fill = "#FFE4E1",
                                                               colour = "#FFE4E1"),
                               legend.background = element_rect(fill = NA),
                               axis.text = element_text(size=15),
                               axis.title = element_text(size=15),
                               title = element_text(size=18),
                               legend.text = element_text(size=18),
                               legend.key = element_rect(color = NA, # filigree color
                                                         fill = NA), # filler color
                               #axis.text.x = element_text(angle=45, hjust=1, vjust=1)
                               legend.key.size = unit(0.3, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 7)

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Supplementary Figure 4 NO.4
#------------------------------------------------------------

fit<-survfit(Surv(OS.time,OS)~Group,data=high)
print(fit)
p4<-ggsurvplot(fit,
               conf.int = FALSE,# Show confidence intervals           
               title = "Stage III+IV", # Setting the title         
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,# Add risk table

               xlab = "Days",# Specify x-axis labels          
               legend = c(0.8,0.9), # Specify legend location     
               legend.title = "", # Setting the legend title       
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_rect(fill = "#FFE4E1",
                                                               colour = "#FFE4E1"),
                               legend.background = element_rect(fill = NA),
                               axis.text = element_text(size=15),
                               axis.title = element_text(size=15),
                               title = element_text(size=18),
                               legend.text = element_text(size=18),
                               legend.key = element_rect(color = NA, # filigree color
                                                         fill = NA), # filler color
                               #axis.text.x = element_text(angle=45, hjust=1, vjust=1)
                               legend.key.size = unit(0.3, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 7)

#######################################   WGD

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Supplementary Figure 4 NO.5
#------------------------------------------------------------

fit<-survfit(Surv(OS.time,OS)~Group,data=WGD)
print(fit)
p5<-ggsurvplot(fit,
               conf.int = FALSE,# Show confidence intervals           
               title = "WGD+", # Setting the title         
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,# Add risk table

               xlab = "Days",# Specify x-axis labels          
               legend = c(0.8,0.9), # Specify legend location     
               legend.title = "", # Setting the legend title       
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_rect(fill = "#FFE4E1",
                                                               colour = "#FFE4E1"),
                               legend.background = element_rect(fill = NA),
                               axis.text = element_text(size=15),
                               axis.title = element_text(size=15),
                               title = element_text(size=18),
                               legend.text = element_text(size=18),
                               legend.key = element_rect(color = NA, # filigree color
                                                         fill = NA), # filler color
                               #axis.text.x = element_text(angle=45, hjust=1, vjust=1)
                               legend.key.size = unit(0.3, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 7)

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Supplementary Figure 4 NO.6
#------------------------------------------------------------

fit<-survfit(Surv(OS.time,OS)~Group,data=nWGD)
print(fit)
p6<-ggsurvplot(fit,
               conf.int = FALSE,# Show confidence intervals           
               title = "WGD-", # Setting the title         
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,# Add risk table

               xlab = "Days",# Specify x-axis labels          
               legend = c(0.8,0.9), # Specify legend location     
               legend.title = "", # Setting the legend title       
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_rect(fill = "#FFE4E1",
                                                               colour = "#FFE4E1"),
                               legend.background = element_rect(fill = NA),
                               axis.text = element_text(size=15),
                               axis.title = element_text(size=15),
                               title = element_text(size=18),
                               legend.text = element_text(size=18),
                               legend.key = element_rect(color = NA, # filigree color
                                                         fill = NA), # filler color
                               #axis.text.x = element_text(angle=45, hjust=1, vjust=1)
                               legend.key.size = unit(0.3, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 7)

#-----------------------------------------------------------
#-----------------------------------------------------------

#######################################   SubType

#------------------------------------------------------------
#    Supplementary Figure 4 NO.7
#------------------------------------------------------------

fit<-survfit(Surv(OS.time,OS)~Group,data=Basal)
print(fit)
p7<-ggsurvplot(fit,
               conf.int = FALSE,# Show confidence intervals           
               title = "Basal", # Setting the title         
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,# Add risk table

               xlab = "Days",# Specify x-axis labels          
               legend = c(0.8,0.9), # Specify legend location     
               legend.title = "", # Setting the legend title       
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_rect(fill = "#FFE4E1",
                                                               colour = "#FFE4E1"),
                               legend.background = element_rect(fill = NA),
                               axis.text = element_text(size=15),
                               axis.title = element_text(size=15),
                               title = element_text(size=18),
                               legend.text = element_text(size=18),
                               legend.key = element_rect(color = NA, # filigree color
                                                         fill = NA), # filler color
                               #axis.text.x = element_text(angle=45, hjust=1, vjust=1)
                               legend.key.size = unit(0.3, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 7)

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Supplementary Figure 4 NO.8
#------------------------------------------------------------

fit<-survfit(Surv(OS.time,OS)~Group,data=Her2)
print(fit)
p8<-ggsurvplot(fit,
               conf.int = FALSE,# Show confidence intervals           
               title = "Her2", # Setting the title         
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,# Add risk table

               xlab = "Days",# Specify x-axis labels          
               legend = c(0.8,0.9), # Specify legend location     
               legend.title = "", # Setting the legend title       
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_rect(fill = "#FFE4E1",
                                                               colour = "#FFE4E1"),
                               legend.background = element_rect(fill = NA),
                               axis.text = element_text(size=15),
                               axis.title = element_text(size=15),
                               title = element_text(size=18),
                               legend.text = element_text(size=18),
                               legend.key = element_rect(color = NA, # filigree color
                                                         fill = NA), # filler color
                               #axis.text.x = element_text(angle=45, hjust=1, vjust=1)
                               legend.key.size = unit(0.3, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 7)

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Supplementary Figure 4 NO.9
#------------------------------------------------------------

fit<-survfit(Surv(OS.time,OS)~Group,data=LumA)
print(fit)
p9<-ggsurvplot(fit,
               conf.int = FALSE,# Show confidence intervals           
               title = "LumA", # Setting the title         
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,# Add risk table

               xlab = "Days",# Specify x-axis labels          
               legend = c(0.8,0.9), # Specify legend location     
               legend.title = "", # Setting the legend title       
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
                               panel.grid.major = element_blank(), 
                               panel.grid.minor = element_blank(),
                               panel.background = element_rect(fill = "#FFE4E1",
                                                               colour = "#FFE4E1"),
                               legend.background = element_rect(fill = NA),
                               axis.text = element_text(size=15),
                               axis.title = element_text(size=15),
                               title = element_text(size=18),
                               legend.text = element_text(size=18),
                               legend.key = element_rect(color = NA, # filigree color
                                                         fill = NA), # filler color
                               #axis.text.x = element_text(angle=45, hjust=1, vjust=1)
                               legend.key.size = unit(0.3, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 7)

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Supplementary Figure 4 NO.10
#------------------------------------------------------------

fit<-survfit(Surv(OS.time,OS)~Group,data=LumB)
print(fit)
p10<-ggsurvplot(fit,
                conf.int = FALSE,# Show confidence intervals           
                title = "LumB", # Setting the title         
                linetype = "strata", # Automatic setting of curve types according to WGD grouping
                palette=c("#EE2C2C","#27408B"), 
                risk.table = TRUE,# Add risk table
 
                xlab = "Days",# Specify x-axis labels          
                legend = c(0.8,0.9), # Specify legend location     
                legend.title = "", # Setting the legend title       
                ggtheme = theme(axis.line =  element_line(colour = "black"),
                                axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
                                panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank(),
                                panel.background = element_rect(fill = "#FFE4E1",
                                                                colour = "#FFE4E1"),
                                legend.background = element_rect(fill = NA),
                                axis.text = element_text(size=15),
                                axis.title = element_text(size=15),
                                title = element_text(size=18),
                                legend.text = element_text(size=18),
                                legend.key = element_rect(color = NA, # filigree color
                                                          fill = NA), # filler color
                                #axis.text.x = element_text(angle=45, hjust=1, vjust=1)
                                legend.key.size = unit(0.3, "inches")
                ), # Setting up the ggplot2 theme
                legend.labs = c("High risk", "Low risk"),
                pval = TRUE,
                pval.size = 7)

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Supplementary Figure 4 NO.11
#------------------------------------------------------------

fit<-survfit(Surv(OS.time,OS)~Group,data=Normal)
print(fit)
p11<-ggsurvplot(fit,
                conf.int = FALSE,# Show confidence intervals           
                title = "Normal", # Setting the title         
                linetype = "strata", # Automatic setting of curve types according to WGD grouping
                palette=c("#EE2C2C","#27408B"), 
                risk.table = TRUE,# Add risk table
 
                xlab = "Days",# Specify x-axis labels          
                legend = c(0.8,0.9), # Specify legend location     
                legend.title = "", # Setting the legend title       
                ggtheme = theme(axis.line =  element_line(colour = "black"),
                                axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
                                panel.grid.major = element_blank(), 
                                panel.grid.minor = element_blank(),
                                panel.background = element_rect(fill = "#FFE4E1",
                                                                colour = "#FFE4E1"),
                                legend.background = element_rect(fill = NA),
                                axis.text = element_text(size=15),
                                axis.title = element_text(size=15),
                                title = element_text(size=18),
                                legend.text = element_text(size=18),
                                legend.key = element_rect(color = NA, # filigree color
                                                          fill = NA), # filler color
                                #axis.text.x = element_text(angle=45, hjust=1, vjust=1)
                                legend.key.size = unit(0.3, "inches")
                ), # Setting up the ggplot2 theme
                legend.labs = c("High risk", "Low risk"),
                pval = TRUE,
                pval.size = 7)


#-----------------------------------------------------------
#-----------------------------------------------------------
