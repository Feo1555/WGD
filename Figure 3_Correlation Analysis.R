setwd("")
library(psych)
library(tidyverse)
library(ggcorrplot)
genes1<-read.table("22beta.txt",header = T,row.names = 1)
genes<-sort(rownames(genes1))

# "brca_247genes_expr" is the expression data of 247 genes in BRCA
brca_247genes_expr<-read.csv("brca_247genes_expr.csv",row.names = 1)
brca_22genes_expr<-brca_247genes_expr[,colnames(brca_247genes_expr) %in% genes]

# "brca_expr" is the expression data of 19718 genes in BRCA
brca_expr<-read.csv(file = "BRCA_expr.csv",row.names = 1)

CCNE2<-brca_expr$CCNE2
KIF18A<-brca_expr$KIF18A
CK<-brca_expr[colnames(brca_expr)%in% c("CCNE2","KIF18A")]
cor.result<-corr.test(brca_22genes_expr,CK,method = "pearson")
cor <- corr.test(CK,brca_22genes_expr,method = "pearson")
cmt <- t(cor$r)
pmt <- t(cor$p)
#------------------------------------------------------------
#    Figure 3A
#------------------------------------------------------------
ggcorrplot(cmt,method = "square",
           outline.color = "white",
           ggtheme = theme_bw(),
           colors = c("#1874CD", "white", "#CC3300"),
           lab = T,
           lab_size=4,
           p.mat=pmt,
           insig="pch",
           pch.col = "red", 
           pch.cex = 6, 
           tl.cex = 9)+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,face = 2),
        axis.text.y = element_text(face = 2),
        axis.ticks= element_blank(),
        panel.border = element_blank(),
        panel.grid.major=element_blank())

#-----------------------------------------------------------
#-----------------------------------------------------------




############################################################

#####KIF18A
KIF18A<-brca_expr$KIF18A
cor_KIF18A<-as.data.frame(matrix(0,22,2))
colnames(cor_KIF18A)<-c("p_value","cor")
rownames(cor_KIF18A)<-colnames(brca_22genes_expr)
for (i in 1:22) {
  a<-cor.test(KIF18A, brca_22genes_expr[,i], method = "pearson")
  cor_KIF18A$p_value[i] <- a$p.value
  cor_KIF18A$cor[i] <- a$estimate
}
#write.csv(cor_KIF18A,file = "cor_KIF18A_22genes.csv",row.names = T)

set.seed(1)
a<-sort(sample(1:401, 281, replace = FALSE))
set.seed(2)
b<-sort(sample(1:586, 410, replace = FALSE))

### training,test,all exp
expr1<-subset(brca_247genes_expr, WGD %in% 1)
expr0<-subset(brca_247genes_expr, WGD %in% 0)
training<-rbind(expr1[a,],expr0[b,])
test<-rbind(expr0[-b,],expr1[-a,])
all<-rbind(expr0,expr1)

expr1<-subset(brca_expr, WGD %in% 1)
expr0<-subset(brca_expr, WGD %in% 0)
training_allgene<-rbind(expr1[a,],expr0[b,])
test_allgene<-rbind(expr0[-b,],expr1[-a,])
all_allgene<-rbind(expr0,expr1)

###all
brca_22genes_expr<-all[,colnames(all) %in% genes]
table(rownames(all) == rownames(brca_22genes_expr))
all$WGD_group[all$WGD %in% 0]<-"WGD-"
all$WGD_group[all$WGD %in% 1]<-"WGD+"
WGD_group<-all[,250]
brca_22genes_expr<-cbind(WGD_group,brca_22genes_expr)
###
exp_rs<-as.data.frame(matrix(NA,987,1))
rownames(exp_rs)<-rownames(all)
colnames(exp_rs)<-"risk_score"
# all risk
for (j in 1:987) {
  a<-0
  for (i in 1:22) {
    b<-as.numeric(genes1[i,1]) * 
      as.numeric(all[j,colnames(all) %in% rownames(genes1)[i]])
    a<-a+b
  }
  exp_rs[j,1]<-a
}
median<-median(exp_rs$risk_score)
exp_rs$risk_group[exp_rs$risk_score <  median] <- "Low risk"
exp_rs$risk_group[exp_rs$risk_score >= median] <- "High risk"
risk_group<-exp_rs$risk_group
table(rownames(exp_rs) == rownames(brca_22genes_expr))
brca_22genes_expr<-cbind(risk_group,brca_22genes_expr)
KIF18A<-all_allgene$KIF18A
brca_23genes_expr<-cbind(KIF18A,brca_22genes_expr)
all_brca_23genes_expr<-brca_23genes_expr

#------------------------------------------------------------
#    Figure 3B NO.1(left to right)
#------------------------------------------------------------

matr<-brca_23genes_expr[,c(1:3,4)]
colnames(matr)[4]<-"exp"
p<-ggplot(matr, aes(KIF18A,exp))+
  geom_point(aes(color=WGD_group))+
  scale_color_manual(values = c("blue","red"),)+
  geom_smooth(method = "lm")+
  stat_cor(data=matr, method = "pearson",size = 15)+
  ylab(colnames(brca_23genes_expr)[4])+
  theme_par()+
  theme(text = element_text(size = 50),legend.position = "bottom",legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(size=10)))

p<-ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
plot_grid(p,ncol = 1,align = "h")

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Figure 3B NO.2
#------------------------------------------------------------

matr<-brca_23genes_expr[,c(1:3,4)]
colnames(matr)[4]<-"exp"
p<-ggplot(matr, aes(KIF18A,exp))+
  geom_point(aes(color=risk_group))+
  scale_color_manual(values = c("red","blue"),)+
  geom_smooth(method = "lm")+
  stat_cor(data=matr, method = "pearson",size = 15)+
  ylab(colnames(brca_23genes_expr)[4])+
  theme_par()+
  theme(text = element_text(size = 50),legend.position = "bottom",legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(size=10)))

p<-ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
plot_grid(p,ncol = 1,align = "h")

#-----------------------------------------------------------
#-----------------------------------------------------------

#####CCNE2
brca_22genes_expr<-brca_247genes_expr[,colnames(brca_247genes_expr) %in% genes]
CCNE2<-brca_expr$CCNE2
cor_CCNE2<-as.data.frame(matrix(0,22,2))
colnames(cor_CCNE2)<-c("p_value","cor")
rownames(cor_CCNE2)<-colnames(brca_22genes_expr)
for (i in 1:22) {
  a<-cor.test(CCNE2, brca_22genes_expr[,i], method = "pearson")
  cor_CCNE2$p_value[i] <- a$p.value
  cor_CCNE2$cor[i] <- a$estimate
}
#write.csv(cor_CCNE2,file = "cor_CCNE2_22genes.csv",row.names = T)

set.seed(1)
a<-sort(sample(1:401, 281, replace = FALSE))
set.seed(2)
b<-sort(sample(1:586, 410, replace = FALSE))

### training,test,all exp
expr1<-subset(brca_247genes_expr, WGD %in% 1)
expr0<-subset(brca_247genes_expr, WGD %in% 0)
training<-rbind(expr1[a,],expr0[b,])
test<-rbind(expr0[-b,],expr1[-a,])
all<-rbind(expr0,expr1)

expr1<-subset(brca_expr, WGD %in% 1)
expr0<-subset(brca_expr, WGD %in% 0)
training_allgene<-rbind(expr1[a,],expr0[b,])
test_allgene<-rbind(expr0[-b,],expr1[-a,])
all_allgene<-rbind(expr0,expr1)




###################################  scatter plot
library(ggplot2)
library(cowplot)
library(ggpubr)
library(ggthemes)
library(ggExtra)


###all
brca_22genes_expr<-all[,colnames(all) %in% genes]
table(rownames(all) == rownames(brca_22genes_expr))
all$WGD_group[all$WGD %in% 0]<-"WGD-"
all$WGD_group[all$WGD %in% 1]<-"WGD+"
WGD_group<-all[,250]
brca_22genes_expr<-cbind(WGD_group,brca_22genes_expr)
###
exp_rs<-as.data.frame(matrix(NA,987,1))
rownames(exp_rs)<-rownames(all)
colnames(exp_rs)<-"risk_score"
# all risk
for (j in 1:987) {
  a<-0
  for (i in 1:22) {
    b<-as.numeric(genes1[i,1]) * 
      as.numeric(all[j,colnames(all) %in% rownames(genes1)[i]])
    a<-a+b
  }
  exp_rs[j,1]<-a
}
median<-median(exp_rs$risk_score)
exp_rs$risk_group[exp_rs$risk_score <  median] <- "Low risk"
exp_rs$risk_group[exp_rs$risk_score >= median] <- "High risk"
risk_group<-exp_rs$risk_group

table(rownames(exp_rs) == rownames(brca_22genes_expr))
brca_22genes_expr<-cbind(risk_group,brca_22genes_expr)

CCNE2<-all_allgene$CCNE2
brca_23genes_expr<-cbind(CCNE2,brca_22genes_expr)

all_brca_23genes_expr<-brca_23genes_expr




###################################  scatter plot
#------------------------------------------------------------
#    Figure 3B NO.3
#------------------------------------------------------------

matr<-brca_23genes_expr[,c(1:3,4)]
colnames(matr)[4]<-"exp"
p<-ggplot(matr, aes(CCNE2,exp))+
  geom_point(aes(color=WGD_group))+
  scale_color_manual(values = c("blue","red"),)+
  geom_smooth(method = "lm")+
  stat_cor(data=matr, method = "pearson",size = 15)+
  ylab(colnames(brca_23genes_expr)[4])+
  theme_par()+
  theme(text = element_text(size = 50),legend.position = "bottom",legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(size=10)))

p<-ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
plot_grid(p,ncol = 1,align = "h")

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Figure 3B NO.4
#------------------------------------------------------------


matr<-brca_23genes_expr[,c(1:3,4)]
colnames(matr)[4]<-"exp"
p<-ggplot(matr, aes(CCNE2,exp))+
  geom_point(aes(color=risk_group))+
  scale_color_manual(values = c("red","blue"),)+
  geom_smooth(method = "lm")+
  stat_cor(data=matr, method = "pearson",size = 15)+
  ylab(colnames(brca_23genes_expr)[i+3])+
  theme_par()+
  theme(text = element_text(size = 50),legend.position = "bottom",legend.title = element_blank())+
  guides(color=guide_legend(override.aes = list(size=10)))

p<-ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE)
plot_grid(p,ncol = 1,align = "h")

#-----------------------------------------------------------
#-----------------------------------------------------------

