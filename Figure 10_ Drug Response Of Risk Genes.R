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

#------------------------------------------------------------
#    Figure 10 A
#------------------------------------------------------------

ggplot()+
  geom_tile(data=new_df1,
            aes(x=Gene,y=Medcine,fill=p_value,alpha=p_value))+
  scale_fill_manual(values = c("white","#c0c0c0",
                               "#808080","#3f3f3f"),
                    label=c(">0.05",
                            "0.01~0.05",
                            "0.001~0.01",
                            "<0.01"))+
  scale_alpha_manual(values = c(0,1,1,1))+
  guides(alpha=F)+
  theme_bw()+
  theme(legend.key = element_rect(colour="black"),
        axis.text.x = element_text(angle = 90,
                                   hjust=1,
                                   vjust=0.5),)+
  coord_equal()+
  geom_point(data=new_df2,
             aes(x=Gene,y=Medcine,
                 size=abs_cor,
                 color=value))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c("#0092FF",'white','red'))

#-----------------------------------------------------------
#-----------------------------------------------------------

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

#------------------------------------------------------------
#    Figure 10B
#------------------------------------------------------------

ggboxplot(auc7,
          x = "Risk Group",
          y = "AUC",
          fill = "Risk Group",
          #notch=TRUE,
          palette = c("#FA7F6F","#82B0D2"))+
  geom_jitter(width = 0.25,
              size = 1,
              #color = "gray",
              alpha=0.5) +
  #geom_point(size = 1)+
  labs(x="")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 4)+ 
  facet_wrap(.~cpd_name,scales = "free",nrow=2)+
  theme(axis.line  = element_blank(),
        strip.background = element_rect(color = "white"),
        panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "#F8F8FF",colour = "#F8F8FF"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        #legend.title=element_blank(),
        legend.position = "top",
        legend.text = element_text(size=13),
        text = element_text(size = 13),
        axis.text = element_text(size=10),
        axis.title = element_text(size=20),
        legend.key.size = unit(0.25, "inches")
  )

#-----------------------------------------------------------
#-----------------------------------------------------------


