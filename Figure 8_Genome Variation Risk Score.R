library(openxlsx)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(ggthemes)
setwd("")

##############################################  train
brca<-read.xlsx("Supplementary File-data1.xlsx",sheet = 1,rowNames = T)
brca<-as.data.frame(t(brca))
rownames(brca)<-brca[,1]
brca<-brca[,-1]
brca_aneuploidy<-brca[,c(1,46)]
brca_aneuploidy$Genome.doublings[brca_aneuploidy$Genome.doublings == 0]<-"nWGD"
brca_aneuploidy$Genome.doublings[brca_aneuploidy$Genome.doublings == 1]<-"WGD"
brca_aneuploidy$array<-substr(rownames(brca_aneuploidy),1,15)
brca_aneuploidy$burden<-as.numeric(brca_aneuploidy$burden)
###  risk score
# "clin_riskscore_train.csv" is organized clinical information on the training cohort
clin_rs_train<-read.csv("clin_riskscore_train.csv",row.names = 1)
clin_rs_train<-clin_rs_train[,c(2,10)]
brca_aneuploidy_rs<-merge(brca_aneuploidy,clin_rs_train,by = "array")
brca_aneuploidy_rs$Group[brca_aneuploidy_rs$group_risk %in% "Low"]<-"Low risk"
brca_aneuploidy_rs$Group[brca_aneuploidy_rs$group_risk %in% "High"]<-"High risk"

#------------------------------------------------------------
#    Figure 8A NO.1 (Left to right)
#------------------------------------------------------------

brca_aneuploidy_rs %>%## 确定x,y
  ggplot(aes(x=group_risk, y=burden,fill=group_risk)) +
  scale_fill_manual(values = c("#56B4E9", "#E69F00"))+ #Sets the color of the fill.
  #Draw a violin map, "color=" set the color of the outline line of the violin map (the following set to the background of white, in fact, that is, do not want to outline the line)
  geom_violin(trim=FALSE,color="white") + 
  #If "trim" is TRUE (the default), trim the end of the violin to the data range. If FALSE, the tail is not trimmed.
  geom_boxplot(width=0.2,position=position_dodge(0.9)) +
  scale_y_continuous(name = "Aneuploidy burden"
                     #,limits = c(-10,30)
  )+
  scale_x_discrete(name = "") +
  theme_par() +
  ggtitle("TCGA-training") +
  theme(text = element_text(size = 50),
        axis.text = element_text(size=50),
        axis.title = element_text(size=50)) +
  guides(fill="none") +
  stat_compare_means(method = "wilcox.test",size = 20,label = "p")#4.8e-10

#-----------------------------------------------------------
#-----------------------------------------------------------

###################################### training exp risk score

beta<-read.csv("beta.txt",row.names = 1,check.name=F)
brca_247genes_expr<-read.csv("brca_247genes_expr.csv",row.names = 1)
expr1<-subset(brca_247genes_expr, WGD %in% 1)
expr0<-subset(brca_247genes_expr, WGD %in% 0)
0.7*nrow(expr1)#280.7
0.7*nrow(expr0)#410.2
set.seed(1)
a<-sort(sample(1:401, 281, replace = FALSE))
set.seed(2)
b<-sort(sample(1:586, 410, replace = FALSE))
### training
training<-rbind(expr0[b,],expr1[a,])
###
exp_rs<-as.data.frame(matrix(NA,691,1))
rownames(exp_rs)<-rownames(training)
colnames(exp_rs)<-"risk_score"
# training risk
for (j in 1:691) {
  a<-0
  for (i in 1:22) {
    b<-as.numeric(beta[i,1]) * 
      as.numeric(training[j,colnames(training) %in% rownames(beta)[i]])
    a<-a+b
  }
  exp_rs[j,1]<-a
}
exp_rs$array<-substr(rownames(exp_rs),1,15)

#######################################################   aneuploidy
aneuploidy<-read.xlsx("Supplementary File-data1.xlsx",sheet = 3,rowNames = T)
colnames(aneuploidy)<-aneuploidy[1,]
aneuploidy<-aneuploidy[-1,]
rownames(aneuploidy)[1]<-"WGD"
aneuploidy<-as.data.frame(t(aneuploidy[c(1,24),]))
aneuploidy$array<-substr(rownames(aneuploidy),1,12)
exp_rs<-as.data.frame(matrix(NA,691,1))
rownames(exp_rs)<-rownames(training)
colnames(exp_rs)<-"risk_score"
# training risk
for (j in 1:691) {
  a<-0
  for (i in 1:22) {
    b<-as.numeric(beta[i,1]) * 
      as.numeric(training[j,colnames(training) %in% rownames(beta)[i]])
    a<-a+b
  }
  exp_rs[j,1]<-a
}

#Take the intersection of the train risk and aneuploidy samples 690
exp_rs$array<-substr(rownames(exp_rs),1,12)
train_exp_rs<-exp_rs[exp_rs$array %in% aneuploidy$array,]
train_aneuploidy<-aneuploidy[aneuploidy$array %in% exp_rs$array,]

#Sorting so that trains correspond to aneuploidy samples
train_exp_rs<-train_exp_rs[order(train_exp_rs$array),]
train_aneuploidy<-train_aneuploidy[order(train_aneuploidy$array),]
train_exp_rs$array == train_aneuploidy$array

############################################################ scatter plot
###
train_aneuploidy_rs<-merge(train_aneuploidy,train_exp_rs,by= "array")

train_aneuploidy_rs$WGD[train_aneuploidy_rs$WGD == 0] <- "WGD-"
train_aneuploidy_rs$WGD[train_aneuploidy_rs$WGD == 1] <- "WGD+"
train_aneuploidy_rs$WGD <- factor(train_aneuploidy_rs$WGD, levels=c('WGD+', 'WGD-'))

train_aneuploidy_rs$burden <- as.numeric(train_aneuploidy_rs$burden)
###
#------------------------------------------------------------
#    Figure 8B NO.1 (Left to right)
#------------------------------------------------------------

p1<-ggplot(train_aneuploidy_rs, aes(risk_score,burden))+
  geom_point(aes(color=WGD))+
  scale_color_manual(values = c("red","blue"),)+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson",
           size = 15,
           cor.coef.name = c("R","Rho","tau"),
           output.type = "latex")+
  ggtitle("TCGA-training")+
  theme_par()+
  theme(text = element_text(size = 50),
        legend.position = "bottom",
        legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(size=10)))+
  ylab("Anueploidy burden")+
  xlab("Risk Score")
ggMarginal(p1, type = "density", groupColour = TRUE, groupFill = TRUE)

#-----------------------------------------------------------
#-----------------------------------------------------------
