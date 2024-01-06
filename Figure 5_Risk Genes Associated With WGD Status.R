library(paletteer)
library(psych)
library(tidyverse)
setwd("")

clin_rs_all<-read.csv("clinic_entire.csv",row.names = 1)

# "brca_247genes_expr" is the expression data of 247 genes in BRCA
brca_247genes_expr<-read.csv("brca_247genes_expr.csv",row.names = 1)
beta<-read.table("22beta.txt")
gene<-rownames(beta)
brca_22<-brca_247genes_expr[,colnames(brca_247genes_expr) %in% gene]
brca_22$array<-substr(rownames(brca_22), 1, 15)

clin_rs_all<-clin_rs_all[which(clin_rs_all$array %in% brca_22$array),]
brca_22<-brca_22[which(brca_22$array %in% clin_rs_all$array),]
brca_22<-brca_22[,-23]
rs<-clin_rs_all$risk_score

cor.result<-corr.test(clin_rs_all$Risk.Score,brca_22,method = "pearson")
cor.result$p %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>% 
  mutate(p_value=case_when(
    value > 0.05 ~ "A",
    value >0.01 & value <= 0.05 ~ "B",
    value > 0.001 & value <= 0.01 ~ "D",
    value <= 0.001 ~ "E")) -> new_df1
colnames(new_df1)<-c("Risk Score","Gene","value","p_value")
new_df1$`Risk Score`<-"Risk Score"

cor.result$r %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>% 
  mutate(abs_cor=abs(value)) -> new_df2
colnames(new_df2)<-c("Risk Score","Gene","value","abs_cor")
new_df2$`Risk Score`<-"Risk Score"

new_df1<-new_df1[22:1,]
new_df2<-new_df2[22:1,]

new_df1$Gene<-factor(new_df1$Gene, levels=new_df1$Gene)
new_df2$Gene<-factor(new_df2$Gene, levels=new_df2$Gene)

#------------------------------------------------------------
#    Figure 5B
#------------------------------------------------------------

ggplot()+
  geom_tile(data=new_df1,
            aes(x=`Risk Score`,y=Gene,fill=p_value,alpha=p_value))+
  scale_fill_manual(values = c("white"),
                    label=c("<0.01"))+
  scale_alpha_manual(values = c(0,1,1,1))+
  guides(alpha=F)+
  theme_bw()+
  theme(legend.key = element_rect(colour="black"),
        axis.text.x = element_text(angle = 90,
                                   hjust=1,
                                   vjust=0.5),)+
  coord_equal()+
  geom_point(data=new_df2,
             aes(x=`Risk Score`,y=Gene,
                 size=abs_cor,
                 color=value))+
  ylab("")+
  xlab("")+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c("#0092FF",'white','red'))

#-----------------------------------------------------------
#-----------------------------------------------------------

############################    Pheatmap
library(pheatmap)
brca_22<-brca_247genes_expr[,colnames(brca_247genes_expr) %in% gene]
brca_22$array<-substr(rownames(brca_22), 1, 15)
risk = clin_rs_all[,c(1,8)]
data<-merge(risk,brca_22)
rownames(data)<-data$array
data<-data[order(data$Risk.Score),]

risk = data[,c(1,2)]
annotation_col<-risk[,2]
annotation_col<-as.data.frame(annotation_col)
rownames(annotation_col)<-risk$array
colnames(annotation_col)<-"risk score"

test<-data[,-c(1,2)]
test2<-as.data.frame(apply(test,2,scale))
rownames(test2)<-rownames(test)
b<-c(1)
c<-c(-1)
test2[test2>b]<-b
test2[test2<c]<-c
test2<-as.data.frame(t(test2))

test2<-test2[,rownames(annotation_col)]
bk <- seq(-1,1,0.01)


#------------------------------------------------------------
#    Figure 5A
#------------------------------------------------------------

p<-pheatmap(test2,
            annotation_col = annotation_col, 
            cluster_cols = F,
            cluster_rows = F,
            #annotation_colors = ann_colors,
            color = c(colorRampPalette(colors = c("#4e8db6","white"))(length(bk)/2),
                      colorRampPalette(colors = c("white","#c65955"))(length(bk)/2)),
            show_colnames = F,
            show_rownames = T,
            treeheight_row = F,
            treeheight_col = F)

#-----------------------------------------------------------
#-----------------------------------------------------------


##########################Data processing before COX
library(survminer)
library(survival)
library(forestplot)
setwd("")

cox_NES<-read.csv("247-122gene-cox.csv",row.names = 1)
sig_cox=cox_NES[cox_NES[,6]<0.05,]#Filtering the training set for p-values less than 0.05

################################################################################ forest map
beta<-read.csv("beta_22genes.csv")
datac <- sig_cox[rownames(sig_cox) %in% beta[,1], ]
genes<-rownames(datac)
count<-rep("N=682", 22)
A<-cbind(genes,count)
datac <- cbind(A,datac)
datac<-as.data.frame(datac)
datac$HR<-as.numeric(datac$HR)
datac$Lower.95<-as.numeric(datac$Lower.95)
datac$Upper.95<-as.numeric(datac$Upper.95)
datac$p_value<-as.numeric(datac$p_value)
datac$p_value<-signif(datac$p_value, digits = 2)

datac$p[datac$p_value <= 0.05 & datac$p_value > 0.01]<-"*"
datac$p[datac$p_value <= 0.01 & datac$p_value > 0.001]<-"**"
datac$p[datac$p_value <= 0.001]<-"***"
datac$p.value<-paste(datac$p_value,datac$p,sep = "  ")
datac$Hazard_ratio<-paste(sprintf("%0.3f", datac$HR), "(", 
                          sprintf("%0.3f", datac$Lower.95), "-", 
                          sprintf("%0.3f", datac$Upper.95), ")", sep = "")

## The rest of the columns in the table.
tabletextc <- cbind(c("Genes",datac$genes),
                    c("Hazard ratio(95% CI)",datac$Hazard_ratio),
                    c("P Value",datac$p.value))


##Forest mapping

#------------------------------------------------------------
#    Figure 5C
#------------------------------------------------------------

forestplot(labeltext=tabletextc,
           graph.pos=2, #Where the Pvalue box line chart is located
           mean=c(NA,datac$HR),
           lower=c(NA,datac$Lower.95), 
           upper=c(NA,datac$Upper.95),
           #Define Title
           title="",
           #Define the x-axis
           #xlab="    <---PCI Better---   ---Medical Therapy Better--->",
           #Setting the line pattern according to the position of the subgroups, the width creates a "blockiness".
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                           "2" = gpar(lty=2),
                           "24"= gpar(lty=2) ),
           #The cex parameter in the fpTxtGp function sets the size of each component
           txt_gp=fpTxtGp(label=gpar(cex=2.0),
                          ticks=gpar(cex=2.0),
                          xlab=gpar(cex = 2.0),
                          title=gpar(cex = 0.7)),
           #The fpColors function sets the colors
           col=fpColors(lines="black", zero = "gray50",box = "red"),
           #Location of baseline in box plot
           zero=1,cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #Box size, line width
           lwd.ci=2, boxsize=0.5,
           #Add small vertical lines at the ends of the box line chart, height
           ci.vertices=F, ci.vertices.height = 0.4)

#-----------------------------------------------------------
#-----------------------------------------------------------

################################################################################ KM
library(survminer)
library(survival)
library(rms)
library(ggrisk)
library(survivalROC)
#################################   training
#"training.RData" is a training group of 22gene expression data and clinical data
load("training.RData")
fit<-survfit(Surv(OS.time,OS)~group,data=dat2)
print(fit)
#------------------------------------------------------------
#    Figure 5D NO.1
#------------------------------------------------------------

p1<-ggsurvplot(fit,
               conf.int = F,# Show confidence intervals          
               title = "Training OS", # Setting the title         
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,# Add risk table
               xlab = "Days",# Specify x-axis labels         
               legend = c(0.8,0.9), # Specify legend location     
               legend.title = "", # Specify legend location   
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
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
                               legend.key.size = unit(0.4, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 9)

#-----------------------------------------------------------
#-----------------------------------------------------------

ddist <- datadist(dat1)
options(datadist="ddist")
fit <- coxph(Surv(OS.time, OS) ~ ANLN+BIRC3+C9orf68+CAMK4+CD200R1+
               CD226+CD79A+CLEC10A+DTHD1+EOMES+
               FCER2+FGL2+KLRB1+LAT+LGALS2+
               NCR3+PPP2R2B+PSTPIP1+PTPRC+SLA2+
               SPN+STAP1,
             data = dat1 )

#------------------------------------------------------------
#    Figure 5D NO.2
#------------------------------------------------------------

#Risk Assessment Scatterplot
p2<-two_scatter(fit, 
                new.data=dat1,#The previous cox regression model
                cutoff.value= 'median' ,#High-risk and low-risk groups by median value
                cutoff.x = 350, #Position of cutoff in Figure A
                cutoff.y = -4,
                code.0 = 'Alive',#The text displayed in Legend in Figure B
                code.1 = 'Dead',
                code.highrisk = 'High' ,#Text displayed by Legend in Figure C
                code.lowrisk = 'Low',
                title.A.ylab='Risk score',#The y-axis labels of Figure A
                title.B.ylab='Survival Time',#The y-axis labels of Figure B
                title.A.legend= 'Group',#Figure A Legend's name
                color.A = c(low = "#27408B", high = "#EE2C2C"),
                color.B = c(code.0 = "#27408B", code.1 = "#EE2C2C"),
                size.AB = 0, #Remove A and B
                size.Atext = 18,
                size.Btext = 18,
                size.ylab.title = 22,
                size.xlab.title = 22,
                size.xtext = 18,
                size.cutoff = 8,
                size.legendtitle = 21,
                size.legendtext = 19
)#Figure C Legend's name

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Figure 5D NO.3
#------------------------------------------------------------

par(mfrow = c(3,1))
nobs <- NROW(dat2)
surv_roc3=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*3,  method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc3$FP, surv_roc3$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="seagreen2" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,cex.main = 1.5,
     xlab="1-Specificity",ylab="Sensitivity",main = "Training OS")

par(new=TRUE)
surv_roc5=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*5, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc5$FP, surv_roc5$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="deepskyblue" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc7=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*7, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc7$FP, surv_roc7$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="yellow" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc10=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*10, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc10$FP, surv_roc10$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="sienna1" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")
lines(x=-1:2,y=-1:2,lwd=2,col="grey60")


auc3=round(surv_roc3$AUC,2)
auc5=round(surv_roc5$AUC,2)
auc7=round(surv_roc7$AUC,2)
auc10=round(surv_roc10$AUC,2)

legend("topleft",inset=c(0.3,0.6),
       c(paste("3-year AUC(",auc3,")",sep = ""),
         paste("5-year AUC(",auc5,")",sep = ""),
         paste("7-year AUC(",auc7,")",sep = ""),
         paste("10-year AUC(",auc10,")",sep = "")), 
       col=c("seagreen2", "deepskyblue","yellow","sienna1"), 
       cex = 1.1,lty =1,lwd = 2, bty = "n")

#-----------------------------------------------------------
#-----------------------------------------------------------

#################################   test
#"test.RData" is a training group of 22gene expression data and clinical data
load("test.RData")
fit<-survfit(Surv(OS.time,OS)~group,data=dat2)
print(fit)

#------------------------------------------------------------
#    Figure 5E NO.1
#------------------------------------------------------------

p4<-ggsurvplot(fit,
               conf.int = F,# Show confidence intervals          
               title = "Test OS", # Setting the title         
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,# Add risk table
               xlab = "Days",# Specify x-axis labels         
               legend = c(0.8,0.9), # Specify legend location     
               legend.title = "", # Specify legend location   
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
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
                               legend.key.size = unit(0.4, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 9)
#-----------------------------------------------------------
#-----------------------------------------------------------

ddist <- datadist(dat1)
options(datadist="ddist")
fit <- coxph( Surv(OS.time, OS) ~ ANLN+BIRC3+C9orf68+CAMK4+CD200R1+
                CD226+CD79A+CLEC10A+DTHD1+EOMES+
                FCER2+FGL2+KLRB1+LAT+LGALS2+
                NCR3+PPP2R2B+PSTPIP1+PTPRC+SLA2+
                SPN+STAP1,
              data = dat1 )

#------------------------------------------------------------
#    Figure 5E NO.2
#------------------------------------------------------------

p5<-two_scatter(fit, 
                new.data=dat1,#The previous cox regression model
                cutoff.value= 'median' ,#High-risk and low-risk groups by median value
                cutoff.x = 145, #Position of cutoff in Figure A
                cutoff.y = -4,
                code.0 = 'Alive',#The text displayed in Legend in Figure B
                code.1 = 'Dead',
                code.highrisk = 'High' ,#Text displayed by Legend in Figure C
                code.lowrisk = 'Low',
                title.A.ylab='Risk score',#The y-axis labels of Figure A
                title.B.ylab='Survival Time',#The y-axis labels of Figure B
                title.A.legend= 'Group',#Figure A Legend's name
                color.A = c(low = "#27408B", high = "#EE2C2C"),
                color.B = c(code.0 = "#27408B", code.1 = "#EE2C2C"),
                size.AB = 0, #Remove A and B
                size.Atext = 18,
                size.Btext = 18,
                size.ylab.title = 22,
                size.xlab.title = 22,
                size.xtext = 18,
                size.cutoff = 8,
                size.legendtitle = 21,
                size.legendtext = 19
)#Figure C Legend's name

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Figure 5E NO.3
#------------------------------------------------------------

nobs <- NROW(dat2)
surv_roc3=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*3,  method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc3$FP, surv_roc3$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="seagreen2" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,cex.main = 1.5,
     xlab="1-Specificity",ylab="Sensitivity",main = "Test OS")

par(new=TRUE)
surv_roc5=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*5, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc5$FP, surv_roc5$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="deepskyblue" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc7=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*7, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc7$FP, surv_roc7$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="yellow" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc10=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*10, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc10$FP, surv_roc10$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="sienna1" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")
lines(x=-1:2,y=-1:2,lwd=2,col="grey60")


auc3=round(surv_roc3$AUC,2)
auc5=round(surv_roc5$AUC,2)
auc7=round(surv_roc7$AUC,2)
auc10=round(surv_roc10$AUC,2)

legend("topleft",inset=c(0.3,0.6),
       c(paste("3-year AUC(",auc3,")",sep = ""),
         paste("5-year AUC(",auc5,")",sep = ""),
         paste("7-year AUC(",auc7,")",sep = ""),
         paste("10-year AUC(",auc10,")",sep = "")), 
       col=c("seagreen2", "deepskyblue","yellow","sienna1"), 
       cex = 1.1,lty =1,lwd = 2, bty = "n")

#-----------------------------------------------------------
#-----------------------------------------------------------
#################################   entire

#"entire.RData" is a training group of 22gene expression data and clinical data
load("entire.RData")

fit<-survfit(Surv(OS.time,OS)~group,data=dat2)
print(fit)

#------------------------------------------------------------
#    Figure 5F NO.1
#------------------------------------------------------------

p7<-ggsurvplot(fit,
               conf.int = F,# Show confidence intervals          
               title = "Entire OS", # Setting the title         
               linetype = "strata", # Automatic setting of curve types according to WGD grouping
               palette=c("#EE2C2C","#27408B"), 
               risk.table = TRUE,# Add risk table
               xlab = "Days",# Specify x-axis labels         
               legend = c(0.8,0.9), # Specify legend location     
               legend.title = "", # Specify legend location   
               ggtheme = theme(axis.line =  element_line(colour = "black"),
                               axis.title.y = element_text(vjust = -6),#Set vjust to a negative value
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
                               legend.key.size = unit(0.4, "inches")
               ), # Setting up the ggplot2 theme
               legend.labs = c("High risk", "Low risk"),
               pval = TRUE,
               pval.size = 9)


#-----------------------------------------------------------
#-----------------------------------------------------------

ddist <- datadist(dat1)
options(datadist="ddist")
fit <- coxph( Surv(OS.time, OS) ~ ANLN+BIRC3+C9orf68+CAMK4+CD200R1+
                CD226+CD79A+CLEC10A+DTHD1+EOMES+
                FCER2+FGL2+KLRB1+LAT+LGALS2+
                NCR3+PPP2R2B+PSTPIP1+PTPRC+SLA2+
                SPN+STAP1,
              data = dat1 )

#------------------------------------------------------------
#    Figure 5E NO.2
#------------------------------------------------------------

p8<-two_scatter(fit, 
                new.data=dat1,#The previous cox regression model
                cutoff.value= 'median' ,#High-risk and low-risk groups by median value
                cutoff.x = 490, #Position of cutoff in Figure A
                cutoff.y = -4,
                code.0 = 'Alive',#The text displayed in Legend in Figure B
                code.1 = 'Dead',
                code.highrisk = 'High' ,#Text displayed by Legend in Figure C
                code.lowrisk = 'Low',
                title.A.ylab='Risk score',#The y-axis labels of Figure A
                title.B.ylab='Survival Time',#The y-axis labels of Figure B
                title.A.legend= 'Group',#Figure A Legend's name
                color.A = c(low = "#27408B", high = "#EE2C2C"),
                color.B = c(code.0 = "#27408B", code.1 = "#EE2C2C"),
                size.AB = 0, #Remove A and B
                size.Atext = 18,
                size.Btext = 18,
                size.ylab.title = 22,
                size.xlab.title = 22,
                size.xtext = 18,
                size.cutoff = 8,
                size.legendtitle = 21,
                size.legendtext = 19
)#Figure C Legend's name

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Figure 5E NO.3
#------------------------------------------------------------

nobs <- NROW(dat2)
surv_roc3=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*3,  method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc3$FP, surv_roc3$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="seagreen2" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,cex.main = 1.5,
     xlab="1-Specificity",ylab="Sensitivity",main = "Entire OS")

par(new=TRUE)
surv_roc5=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*5, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc5$FP, surv_roc5$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="deepskyblue" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc7=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*7, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc7$FP, surv_roc7$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="yellow" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc10=survivalROC(dat2$OS.time, dat2$OS, dat2$risk_score, predict.time=365*10, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc10$FP, surv_roc10$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="sienna1" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")
lines(x=-1:2,y=-1:2,lwd=2,col="grey60")


auc3=round(surv_roc3$AUC,2)
auc5=round(surv_roc5$AUC,2)
auc7=round(surv_roc7$AUC,2)
auc10=round(surv_roc10$AUC,2)

legend("topleft",inset=c(0.3,0.6),
       c(paste("3-year AUC(",auc3,")",sep = ""),
         paste("5-year AUC(",auc5,")",sep = ""),
         paste("7-year AUC(",auc7,")",sep = ""),
         paste("10-year AUC(",auc10,")",sep = "")), 
       col=c("seagreen2", "deepskyblue","yellow","sienna1"), 
       cex = 1.1,lty =1,lwd = 2, bty = "n")

#-----------------------------------------------------------
#-----------------------------------------------------------
