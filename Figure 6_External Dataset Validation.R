library(stringr)
library(survival)
library(survminer)
library(survivalROC)
setwd("")

####################GSE37751 and GPL6244 are from the GEO database 

exp<-read.csv("GSE37751_series_matrix.txt",skip = 75,sep = "\t")
gpl<-read.csv("GPL6244-17930.txt",skip = 12,sep = "\t",check.names = F)
clinic<-read.csv("GSE37751_clinicalData.txt",sep = "",check.names = F,header = F)
clinic<-clinic[c(2,19,18),-1]
rownames(clinic)<-c("Sample","OS","OS.time")
clinic<-as.data.frame(t(clinic))
clinic$OS<-gsub("survival event: ","",clinic$OS)#(1 = dead, 0 = alive or censored)
clinic$OS[clinic$OS %in% "not applicable"]<-NA
clinic$OS[clinic$OS %in% "Dead"] <- 1
clinic$OS[clinic$OS %in% "Alive"] <- 0
clinic$OS.time<-substr(clinic$OS.time,20,21)
clinic$OS.time[clinic$OS.time %in% "no"]<-NA#months
rownames(clinic)<-1:nrow(clinic)
gpl<-gpl[,c(1,10)]
gpl[,2]<-str_extract(str_extract(gpl[,2],"\\//.+?\\//"),"\\ .+?\\ ") 
gpl[,2]<-gsub(" ","",gpl[,2])
exp1<-merge(gpl,exp,by.x = "ID", by.y = "ID_REF" )
beta<-read.csv("22beta.txt",sep = "")
genes<-rownames(beta)
colnames(exp1)[2]<- "Gene symbol"
exp2<-exp1[exp1$`Gene symbol` %in% genes,]

#Identical genes take the mean of the expression values
a<-as.data.frame(table(exp2$`Gene symbol`))
a1<-subset(a,Freq > 1)
exp2_1<-exp2[!exp2$`Gene symbol` %in% a1$Var1,]
exp2_1<-exp2_1[,-1]
exp2_2<-data.frame()
for (i in 1:nrow(a1)) {
  matr<-exp2[exp2$`Gene symbol` %in% a1$Var1[i],]
  matr1<-as.data.frame(t(colMeans(matr[,-c(1,2)])))
  matr1<-cbind(matr[1,2],matr1)
  exp2_2<-rbind(exp2_2,matr1)
}
colnames(exp2_2)[1]<-colnames(exp2_1)[1]
exp3<-rbind(exp2_1,exp2_2)
rownames(exp3)<-1:nrow(exp3)


#################################   Calculating the risk score

exp4<-as.data.frame(t(exp3))
colnames(exp4)<-exp3[,1]
exp4<-exp4[-1,]
risk_score<-as.data.frame(matrix(NA,108,1))
rownames(risk_score)<-rownames(exp4)
colnames(risk_score)<-"risk_score"

for (j in 1:108) {
  a<-0
  for (i in 1:21) {
    b<-as.numeric(beta[rownames(beta) %in% colnames(exp4)[i],1]) * as.numeric(exp4[j,colnames(exp4)[i]])
    a<-a+b
  }
  risk_score[j,1]<-a
}

rownames(clinic)<-clinic$Sample
rownames(clinic) == rownames(risk_score)
clinic<-cbind(clinic,risk_score)
clinic$OS.time
clinic$OS[clinic$OS %in% "NA"] <- NA
clinic$OS.time[clinic$OS.time %in% "NA"] <- NA
clinic$OS<-as.numeric(clinic$OS)
clinic$OS.time<-as.numeric(clinic$OS.time)
clinic$OS.time<-clinic$OS.time*30

############################################################  KM curve diagram

clinic<-na.omit(clinic)
median<-median(clinic$risk_score)
clinic[clinic$risk_score < median(clinic$risk_score),5]  <- "Low"
clinic[clinic$risk_score >= median(clinic$risk_score),5] <- "High"
colnames(clinic)[5] <- "group_risk"
fit<-survfit(Surv(OS.time,OS)~group_risk,data=clinic)
print(fit)

#------------------------------------------------------------
#    Figure 6 NO.1 (Left to right)
#------------------------------------------------------------

ggsurvplot(fit,
           conf.int = F,# Show confidence intervals          
           title = "GSE37751 OS", # Setting the title         
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
                           axis.text = element_text(size=13),
                           axis.title = element_text(size=13),
                           title = element_text(size=16),
                           legend.text = element_text(size=16),
                           legend.key = element_rect(color = NA, # filigree color
                                                     fill = NA), # filler color
                           legend.key.size = unit(0.4, "inches")
           ), # Setting up the ggplot2 theme
           legend.labs = c("High risk", "Low risk"),
           pval = TRUE,
           pval.size = 7)

#-----------------------------------------------------------
#-----------------------------------------------------------

############################################################  ROC

#------------------------------------------------------------
#    Figure 6 NO.5 (Left to right)
#------------------------------------------------------------

nobs <- NROW(clinic)
surv_roc3=survivalROC(clinic$OS.time, clinic$OS, clinic$risk_score, predict.time=365*3,  method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc3$FP, surv_roc3$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="seagreen2" , lwd=2,
     cex.axis  = 1.5,cex.lab = 1.5,cex.main = 1.5,
     xlab="1-Specificity",ylab="Sensitivity",main = "GSE37751 OS")

par(new=TRUE)
surv_roc5=survivalROC(clinic$OS.time, clinic$OS, clinic$risk_score, predict.time=365*5, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc5$FP, surv_roc5$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="deepskyblue" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc7=survivalROC(clinic$OS.time, clinic$OS, clinic$risk_score, predict.time=365*7, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc7$FP, surv_roc7$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="yellow" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc10=survivalROC(clinic$OS.time, clinic$OS, clinic$risk_score, predict.time=365*10, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc10$FP, surv_roc10$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="sienna1" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")
lines(x=-1:2,y=-1:2,lwd=2,col="grey60")


auc3=round(surv_roc3$AUC,2)
auc5=round(surv_roc5$AUC,2)
auc7=round(surv_roc7$AUC,2)
auc10=round(surv_roc10$AUC,2)

legend("topleft",inset=c(0.28,0.6),
       c(paste("3-year AUC(",auc3,")",sep = ""),
         paste("5-year AUC(",auc5,")",sep = ""),
         paste("7-year AUC(",auc7,")",sep = ""),
         paste("10-year AUC(",auc10,")",sep = "")),
       col=c("seagreen2", "deepskyblue","yellow","sienna1"), 
       cex = 1.3,lty =1,lwd = 2, bty = "n")

#-----------------------------------------------------------
#-----------------------------------------------------------

