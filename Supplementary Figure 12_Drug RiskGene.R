library(psych)
library(tidyverse)
library(paletteer)
library(ggpubr)
library(Hmisc)
setwd("")

# "cl_gene_exp.csv" is integrated to obtain gene expression data for 40 cell lines in rows 
# and 225 genes in columns
cl_gene_exp<-read.csv("cl_gene_exp.csv",row.names = 1,check.name=F)

# "cl_medicine_auc.csv"  is obtained by integrating drug sensitivity data (AUC) 
# for 40 cell lines in rows and 545 drugs in columns.
cl_medicine_auc<-read.csv("cl_medicine_auc.csv",row.names = 1,check.name=F)
beta<-read.csv("beta_22genes.csv",row.names = 1,check.name=F)
gene<-rownames(beta)

# "40medcine.csv"is significantly associated with risk genes for 40 different anticancer drugs
medcine<-read.csv("40medcine.csv",row.names = 1)
medcine<-medcine$x
cl_gene_exp1<-cl_gene_exp[,colnames(cl_gene_exp) %in% gene]
cl_medicine_auc1<-cl_medicine_auc[,colnames(cl_medicine_auc) %in% medcine]
cor.result<-corr.test(cl_gene_exp1,cl_medicine_auc1,method = "pearson")

cor.result$p %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>% 
  mutate(p_value=case_when(
    value > 0.05 ~ "A",
    value >0.01 & value <= 0.05 ~ "B",
    value > 0.001 & value <= 0.01 ~ "D",
    value <= 0.001 ~ "E")) -> new_df1
colnames(new_df1)<-c("Gene","Medcine","value","p_value")
new_df1[which(new_df1$value %in% NA),3]<-0
new_df1[which(new_df1$p_value %in% NA),4]<-"A"

cor.result$r %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>% 
  mutate(abs_cor=abs(value)) -> new_df2
colnames(new_df2)<-c("Gene","Medcine","value","abs_cor")

### difference in AUC values of 7 drug component

# "auc5.csv" is based on paired drug susceptibility data (AUC) obtained 
# by integrating 40 cell lines and 545 drugs. 
auc5<-read.csv("auc5.csv")
auc6<-auc5[,c(3,8,10,27)]
test3<-read.csv("test3.csv")
test4<-test3[,c(1,3)]
test4$risk_group <- test4$risk_score<median(test4$risk_score)
test4$risk_group[test4$risk_group == "TRUE"] <- "Low"
test4$risk_group[test4$risk_group == "FALSE"] <- "High"
auc6<-merge(auc6,test4,by = "DepMap_ID")
colnames(auc6)<-c("DepMap_ID","WGD","cpd_name","AUC","risk_score","risk_group")
auc6<-auc6[order(auc6$cpd_name),]

#######################################       7 drugs with a significant positive correlation
auc7<-auc6[auc6$cpd_name %in% c("AZD7762",
                                "CIL41",
                                "dasatinib",
                                "GSK-3 inhibitor IX",
                                "semagacestat",
                                "staurosporine",
                                "vandetanib"),]
colnames(auc7)[6]<-"Risk Group"
auc7$cpd_name<-capitalize(auc7$cpd_name)

auc7$WGD[auc7$WGD %in% 0]<-"WGD-"
auc7$WGD[!auc7$WGD %in% "WGD-"]<-"WGD+"
auc7$WGD <- factor(auc7$WGD, levels=c('WGD+', 'WGD-'))
#relation

#------------------------------------------------------------
#    Supplementary Figure 12 
#------------------------------------------------------------

ggplot(auc7, aes(risk_score,AUC))+
  geom_point(size=2,aes(color=WGD))+
  geom_smooth(method = "lm")+
  labs(x="Risk Score")+
  stat_cor(method = "pearson",
           label.sep = "\n",
           size = 4)+ 
  facet_wrap(.~cpd_name,scales = "free",nrow=1)+
  scale_discrete_manual(values=c("#FA7F6F","#82B0D2"),
                        aesthetics = 'colour')+
  theme(axis.line  = element_blank(),
        strip.background = element_rect(color = "white"),
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "#F8F8FF",colour = "#F8F8FF"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.title=element_blank(),
        legend.position = "top",
        legend.text = element_text(size=16),
        text = element_text(size = 16),
        axis.text = element_text(size=10),
        axis.title = element_text(size=20),
        legend.key.size = unit(0.25, "inches")
  )

#-----------------------------------------------------------
#-----------------------------------------------------------
