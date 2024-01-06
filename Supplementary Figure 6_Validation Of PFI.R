########################## Data processing before cox
library(survminer)
library(survival)
library(rms)
library(ggrisk)
library(survivalROC)
setwd("")

#"brca_247genes_expr.csv" is the expression data of 247 genes in BRCA from the TCGA database
brca_247genes_expr<-read.csv("brca_247genes_expr.csv",row.names = 1)
expr1<-subset(brca_247genes_expr, WGD %in% 1)
expr0<-subset(brca_247genes_expr, WGD %in% 0)
0.7*nrow(expr1)#280.7
0.7*nrow(expr0)#410.2
set.seed(1)
a<-sort(sample(1:401, 281, replace = FALSE))
set.seed(2)
b<-sort(sample(1:586, 410, replace = FALSE))
training<-rbind(expr0[b,],expr1[a,])

#############################################
#"TCGA_CDR_clinical_data.txt" is clinical data on breast cancer patients from the TCGA database
clin<-read.table("TCGA_CDR_clinical_data.txt",header = T,sep = "\t")
pfi<-clin[,colnames(clin) %in% "PFI" | colnames(clin) %in% "PFI.time"]
rownames(pfi)=clin$bcr_patient_barcod
pfi=na.omit(pfi[-which(pfi$PFI.time==0),])
expr<-training[,-c(1,2)]
name1=substr(rownames(expr),1,12)
rownames(expr)=name1
pfi1=pfi[intersect(name1,rownames(pfi)),]
expr=expr[intersect(name1,rownames(pfi)),]
dat=cbind(expr,pfi1)

#############################COX
#"247-122gene-cox.csv" is the result of cox regression analysis
cox_NES<-read.csv("247-122gene-cox.csv",row.names = 1)
sig_glmname<-c("ANLN","BIRC3","C9orf68","CAMK4","CD200R1",
               "CD226","CD79A","CLEC10A","DTHD1","EOMES",
               "FCER2","FGL2","KLRB1","LAT","LGALS2",
               "NCR3","PPP2R2B","PSTPIP1","PTPRC","SLA2","SPN","STAP1"  )
sig_cox=cox_NES[rownames(cox_NES) %in% sig_glmname,]
dat1=dat[,c(which(colnames(dat) %in% sig_glmname),248,249)]

#################################risk score
sig_training_data=as.matrix(dat[,rownames(sig_cox)])
beta<-read.table("22beta.txt")
risk_score<-as.data.frame(matrix(NA,682,1))
rownames(risk_score)<-rownames(sig_training_data)
colnames(risk_score)<-"risk_score"
for (j in 1:682) {
  a<-0
  for (i in 1:22) {
    b<-as.numeric(beta[rownames(beta) %in% sig_glmname[i],1]) * sig_training_data[j,colnames(sig_training_data) %in% sig_glmname[i]]
    a<-a+b
  }
  risk_score[j,1]<-a
}
dat1=dat[,c(which(colnames(dat) %in% sig_glmname),248,249)]
risk_score$group<-0
risk_score[risk_score$risk_score < median(risk_score$risk_score),2] <- "Low"
risk_score[risk_score$risk_score >= median(risk_score$risk_score),2] <- "High"
risk_score$patient_barcode<-substr(rownames(risk_score),1,12)
dat1$patient_barcode<-rownames(dat1)
dat2<-merge(dat1,risk_score,by  = "patient_barcode")


########################## Data processing before cox
pfi<-clin[,colnames(clin) %in% "PFI" | 
            colnames(clin) %in% "PFI.time"| 
            colnames(clin) %in% "ajcc_pathologic_tumor_stage"| 
            colnames(clin) %in% "age_at_initial_pathologic_diagnosis" ]

rownames(pfi)=clin$bcr_patient_barcod
pfi=na.omit(pfi[-which(pfi$PFI.time==0),])

wgd<-training[,c(1,2)]
name1=substr(rownames(wgd),1,12)
rownames(wgd)=name1
pfi1=pfi[intersect(name1,rownames(pfi)),]
wgd=wgd[intersect(name1,rownames(pfi)),]
data=cbind(wgd,pfi1)
data$patient_barcode <- substr(data$array,1,12)
data1<-merge(data,risk_score,by = "patient_barcode",)
colnames(data1)<-c("patient_barcode","array","WGD",
                   "age","tumor_stage","PFI",
                   "PFI.time","risk_score","group")

table(data1$tumor_stage,data1$PFI)
data1$tumor_stage[data1$tumor_stage %in%  "Stage I" | 
                    data1$tumor_stage %in%  "Stage IA"| 
                    data1$tumor_stage %in%  "Stage IB"] <- 1

data1$tumor_stage[data1$tumor_stage %in%  "Stage II" | 
                    data1$tumor_stage %in%  "Stage IIA"| 
                    data1$tumor_stage %in%  "Stage IIB"] <- 2

data1$tumor_stage[data1$tumor_stage %in%  "Stage III" | 
                    data1$tumor_stage %in%  "Stage IIIA"| 
                    data1$tumor_stage %in%  "Stage IIIB"| 
                    data1$tumor_stage %in%  "Stage IIIC"] <- 3


data1$tumor_stage[data1$tumor_stage %in%  "Stage IV" ]<- 4

data1$tumor_stage[data1$tumor_stage %in%  "Stage X" | 
                    data1$tumor_stage %in%  "[Discrepancy]"| 
                    data1$tumor_stage %in%  "[Not Available]"]<- NA


data1$tumor_stage<-as.integer(data1$tumor_stage)

########################
clin_rs_training<-read.csv("clinic_train.csv",row.names = 1)
clin_rs_training<-clin_rs_training[,c(1,10)]
data3<-merge(data1,clin_rs_training,all.x = T)
colnames(data3)<-c("array","patient_barcode","WGD","Age","Stage",
                   "PFI","PFI.time","Risk Score","Group","SubType")

#载入并查看数据集
data2<-data3
head(data3)
ddist <- datadist(data3)
options(datadist="ddist")
data3$WGD<-as.numeric(data3$WGD)

data3$WGD<-as.numeric(data3$WGD)
data3$Stage<-as.numeric(data3$Stage)
data3$SubType[data3$SubType %in% "Basal"]<-1
data3$SubType[data3$SubType %in% "Her2"]<-2
data3$SubType[data3$SubType %in% "LumA"]<-3
data3$SubType[data3$SubType %in% "LumB"]<-4
data3$SubType[data3$SubType %in% "Normal"]<-5
data3$SubType<-as.numeric(data3$SubType)

#################################   training
fit<-survfit(Surv(PFI.time,PFI)~group,data=dat2)
print(fit)
#------------------------------------------------------------
#    Supplementary Figure 6A NO.1
#------------------------------------------------------------

p1<-ggsurvplot(fit,
               conf.int = F,# Show confidence intervals           
               title = "Training PFI", # Setting the title         
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

####################################################### Risk Assessment Scatterplot
ddist <- datadist(dat)
options(datadist="ddist")
fit <- coxph(Surv(PFI.time, PFI) ~ ANLN+BIRC3+C9orf68+CAMK4+CD200R1+
               CD226+CD79A+CLEC10A+DTHD1+EOMES+
               FCER2+FGL2+KLRB1+LAT+LGALS2+
               NCR3+PPP2R2B+PSTPIP1+PTPRC+SLA2+
               SPN+STAP1,
             dat)

#------------------------------------------------------------
#    Supplementary Figure6A NO.2
#------------------------------------------------------------

p2<-two_scatter(fit, 
                new.data=dat1,#The previous cox regression model
                cutoff.value= 'median' ,#High-risk and low-risk groups by median value
                cutoff.x = 350, #cLocation of utoff in Figure A
                cutoff.y = -2.5,
                code.0 = 'Alive',#The text displayed in Legend in Figure B
                code.1 = 'Dead',
                code.highrisk = 'High' ,# Text displayed by Legend in Figure C
                code.lowrisk = 'Low',
                title.A.ylab='Risk score',#Labeling of the y-axis of #Figure A
                title.B.ylab='Survival Time',#The y-axis labels of Fig. B.
                title.A.legend= 'Group',#Figure A Legend's name
                color.A = c(low = "#27408B", high = "#EE2C2C"),
                color.B = c(code.0 = "#27408B", code.1 = "#EE2C2C"),
                size.AB = 0, #去掉AB
                size.Atext = 18,
                size.Btext = 18,
                size.ylab.title = 22,
                size.xlab.title = 22,
                size.xtext = 18,
                size.cutoff = 8,
                size.legendtitle = 21,
                size.legendtext = 19
                #title.C.legend='Expression'
)#Figure C Legend's name

#-----------------------------------------------------------
#-----------------------------------------------------------



#################################   ROC

#------------------------------------------------------------
#    Supplementary Figure 6A NO.3
#------------------------------------------------------------

nobs <- NROW(dat2)
surv_roc3=survivalROC(dat2$PFI.time, dat2$PFI, dat2$risk_score, predict.time=365*3,  method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc3$FP, surv_roc3$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="seagreen2" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,cex.main = 1.5,
     xlab="1-Specificity",ylab="Sensitivity",main = "Training PFI")

par(new=TRUE)
surv_roc5=survivalROC(dat2$PFI.time, dat2$PFI, dat2$risk_score, predict.time=365*5, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc5$FP, surv_roc5$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="deepskyblue" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,cex.main = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc7=survivalROC(dat2$PFI.time, dat2$PFI, dat2$risk_score, predict.time=365*7, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc7$FP, surv_roc7$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="yellow" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,cex.main = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc10=survivalROC(dat2$PFI.time, dat2$PFI, dat2$risk_score, predict.time=365*10, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc10$FP, surv_roc10$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="sienna1" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,cex.main = 1.5,
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

