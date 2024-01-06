library(openxlsx)
library(data.table)
library(ggplot2)
library(RColorBrewer)
library(plyr)
setwd("")
##############################################  train
kaf<-as.data.frame(matrix(0,3,2))
rownames(kaf)<-c("Training","Test","All")
colnames(kaf)<-c("Statistic","P.Value")

###  risk score
clin_rs_train<-read.csv("clin_riskscore_train.csv",row.names = 1)
risk_WGD<-clin_rs_train[,c(2,3,10)]
setDT(risk_WGD)
risk_WGD[,WGD:=ifelse(risk_WGD$WGD %in% 1,"WGD","nWGD")]
ka<-xtabs(~risk_WGD$WGD+risk_WGD$group_risk,data=risk_WGD)
katrain<-chisq.test(ka)
kaf$Statistic[1]<-katrain$statistic
kaf$P.Value[1]<-katrain$p.value
data<-risk_WGD
data$WGD<-as.factor(data$WGD)
data[,group_risk:=ifelse(group_risk %in% "High","High Risk","Low Risk")]
pvalue <- katrain$p.value #Chi-Squared Test

a<- data.frame(table(data$group_risk,data$WGD))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")

#------------------------------------------------------------
#    Supplementary Figure 10A NO.1
#------------------------------------------------------------

ggplot(a,aes(Var1,percent,fill=Var2))+
  geom_bar(stat="identity",position = position_stack(),width =0.7)+
  scale_fill_manual(values = c("blue","red"),label=c("WGD-","WGD+"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     breaks=seq(0, 100, 20),
                     expand=c(0,0), 
                     limits=c(0, 120))+ #Percentage y-axis
  ggtitle("TCGA-training")+
  labs(x="",y="Sample Percent(%)",fill="")+
  geom_text(aes(label=label),vjust=c(-2.5,4,0.5,3.5),size=20,color="white")+
  annotate(geom = "text",
           cex=15,
           x=1.5, y=110, # Adjust the position of the p value according to your own data
           label=paste0("P ", 
                        ifelse(pvalue<0.001, "< 0.001", 
                               paste0("= ",
                                      round(pvalue,3)))),# Add P-value
           color="black") +  
  theme_classic()+
  theme(legend.position = "right",
        legend.text = element_text(size=50),
        text = element_text(size = 50),
        axis.text = element_text(size=50),
        axis.title = element_text(size=50),
        legend.key.size = unit(1, "inches"))

#-----------------------------------------------------------
#-----------------------------------------------------------

#"mutation_burden_wgd01.xlsx" is the mutation burden for each sample form TCGA
mutation_burden<-read.xlsx("mutation_burden_wgd01.xlsx",sheet = "BRCA")
mutation_burden$array<-substr(mutation_burden$Tumor_Sample_Barcode, 1, 15)
mutation_burden<-mutation_burden[,c(4,2)]
kaf<-as.data.frame(matrix(0,3,2))
rownames(kaf)<-c("Training","Test","All")
colnames(kaf)<-c("Statistic","P.Value")

##############################################  train
###  risk score
clin_rs_train<-read.csv("clin_riskscore_train.csv",row.names = 1)
clin_rs_train<-clin_rs_train[,c(2,3,10)]
brca_mutation_rs<-merge(mutation_burden,clin_rs_train,by = "array",all.y = T)

setDT(brca_mutation_rs)
brca_mutation_rs[,mutation:=ifelse(mutation_burden %in% NA,0,1)]
ka<-xtabs(~brca_mutation_rs$mutation+brca_mutation_rs$group_risk,data=brca_mutation_rs)
katrain<-chisq.test(ka)
kaf$Statistic[1]<-katrain$statistic
kaf$P.Value[1]<-katrain$p.value
kad<-as.data.frame(ka)

data<-brca_mutation_rs
data[,group_risk:=ifelse(group_risk %in% "High","High Risk","Low Risk")]
data[,mutation:=ifelse(mutation %in% "1","Mutation","nMutation")]

a<- data.frame(table(data$group_risk,data$mutation))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
pvalue <- katrain$p.value #Chi-Squared Test

#------------------------------------------------------------
#    Supplementary Figure 10B NO.1
#------------------------------------------------------------

ggplot(a,aes(Var1,percent,fill=Var2))+
  geom_bar(stat="identity",position = position_stack(),width =0.7)+
  scale_fill_manual(values = c("#D9761B","#3574B6"),label=c("Mut+","Mut-"))+
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     breaks=seq(0, 100, 20),
                     expand=c(0,0), 
                     limits=c(0, 120))+ #Percentage y-axis
  ggtitle("TCGA-training")+
  labs(x="",y="Sample Percent(%)",fill="")+
  geom_text(aes(label=label),vjust=c(0,3,3,2),size=20,color="white")+
  annotate(geom = "text",
           cex=15,
           x=1.5, y=110, # Adjust the position of the p value according to your own data
           label=paste0("P ", 
                        ifelse(pvalue<0.001, "< 0.001", 
                               paste0("= ",
                                      round(pvalue,3)))),# add p-value
           color="black") +  
  theme_classic()+
  theme(legend.position = "right",
        legend.text = element_text(size=50),
        text = element_text(size = 50),
        axis.text = element_text(size=50),
        axis.title = element_text(size=50),
        legend.key.size = unit(1, "inches"))

#-----------------------------------------------------------
#-----------------------------------------------------------