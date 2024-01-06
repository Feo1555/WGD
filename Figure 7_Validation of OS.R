library(survminer)
library(survival)
library(rms)
library(tableone)
library(ggDCA)
library(foreign)
library(ggthemes)
setwd("")
########################## Data processing before cox

brca_247genes_expr<-read.csv(file = "brca_247genes_expr.csv",row.names = 1)
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
clin<-read.table("TCGA_CDR_clinical_data.txt",header = T,sep = "\t")
os<-clin[,colnames(clin) %in% "OS" | colnames(clin) %in% "OS.time"]
rownames(os)=clin$bcr_patient_barcod
os=na.omit(os[-which(os$OS.time==0),])
expr<-training[,-c(1,2)]
name1=substr(rownames(expr),1,12)
rownames(expr)=name1
os1=os[intersect(name1,rownames(os)),]
expr=expr[intersect(name1,rownames(os)),]
dat=cbind(expr,os1)

#############################COX
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
os<-clin[,colnames(clin) %in% "OS" | 
           colnames(clin) %in% "OS.time"| 
           colnames(clin) %in% "ajcc_pathologic_tumor_stage"| 
           colnames(clin) %in% "age_at_initial_pathologic_diagnosis" ]

rownames(os)=clin$bcr_patient_barcod
os=na.omit(os[-which(os$OS.time==0),])
wgd<-training[,c(1,2)]
name1=substr(rownames(wgd),1,12)
rownames(wgd)=name1
os1=os[intersect(name1,rownames(os)),]
wgd=wgd[intersect(name1,rownames(os)),]
data=cbind(wgd,os1)
data$patient_barcode <- substr(data$array,1,12)
data1<-merge(data,risk_score,by = "patient_barcode",)
colnames(data1)<-c("patient_barcode","array","WGD",
                   "age","tumor_stage","OS",
                   "OS.time","risk_score","group")
table(data1$tumor_stage,data1$OS)
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

data1<-na.omit(data1)
data1$tumor_stage<-as.integer(data1$tumor_stage)

########################
clin_rs_training<-read.csv("clinic_train.csv",row.names = 1)
clin_rs_training<-clin_rs_training[,c(1,10)]
data3<-merge(data1,clin_rs_training,all.x = T)
colnames(data3)<-c("array","patient_barcode","WGD","Age","Stage",
                   "OS","OS.time","Risk Score","Group","SubType")


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


#cox multifactor regression analysis
mul_cox <- coxph(Surv(OS.time, OS) ~ Age+`Risk Score`+Stage+WGD+SubType, data = data3)
#Extraction: variable + HR + 95% CI + 95% CI
mul_cox1 <- summary(mul_cox)
colnames(mul_cox1$conf.int)
multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 3))
#Extraction: HR (95% CI) and P
multi2<-ShowRegTable(mul_cox, 
                     exp=TRUE, 
                     digits=3, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)
#Combine the results of the two extractions into a table; name the "result"
result <-cbind(multi1,multi2)
#The row name is converted to the first column of the table and given the name "Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result

datam<-result[,c(1,5,6)]
datam[2,1]<-"Risk Score"
rownames(datam)<-NULL
colnames(datam)<-NULL
datam<-as.data.frame(datam)

datam$p<-""
datam$p[1:3]<-"***"
datam$p[4:5]<-""
datam[,3]<-paste(datam[,3],datam$p,sep = "  ")
datam<-datam[,-4]
datam<-rbind(c("","Hazard ratio(95% CI)","P Value"),datam)
colnames(datam)<-NULL
datam<-as.matrix(datam)

datam[,2] <- gsub("[","(",datam[,2],fixed = TRUE)
datam[,2] <- gsub("]",")",datam[,2],fixed = TRUE)
datam[,2] <- gsub(", ","-",datam[,2],fixed = TRUE)

#------------------------------------------------------------
#    Figure 7B
#------------------------------------------------------------

forestplot(labeltext=datam,
           graph.pos=2, #the location of the Pvalue box plot.
           mean=c(NA,result$`exp(coef)`),
           lower=c(NA,result$`lower .95`), 
           upper=c(NA,result$`upper .95`),
           #Define Title
           title="",
           ##Define the x-axis
           #xlab="    <---PCI Better---   ---Medical Therapy Better--->",
           ##Setting the line pattern according to the position of the subgroups, the width creates a "blockiness".
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                           "2" = gpar(lty=2),
                           #"4" = gpar(lwd=30, lineend="butt", columns=c(1:4), col="#99999922"),
                           "7"= gpar(lty=2) ),
           #The cex parameter in the fpTxtGp function sets the size of each component
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           ##The fpColors function sets the colors
           col=fpColors(box = "red", lines="black", zero = "gray50"),
           #Location of baseline in box plot
           zero=1,
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #Box size, line width
           lwd.ci=1, boxsize=0.3,
           #Add small vertical lines at the ends of the box-and-line diagram, height
           ci.vertices=F, ci.vertices.height = 0.4)

#-----------------------------------------------------------
#-----------------------------------------------------------

############################################################ forest map
data2<-data3[,c(4,8,5,3,10,6,7)]
cox_NES<-matrix(nrow=ncol(data2)-2,ncol=6)
for (j in 1:(ncol(data2)-2)) {
  Bcox<-coxph(Surv(OS.time, OS)~data2[,j],data=data2)
  summcph<-summary(Bcox)
  cox_NES[j,1]<-summcph$conf.int[1]
  cox_NES[j,2]<-summcph$conf.int[3]
  cox_NES[j,3]<-summcph$conf.int[4]
  cox_NES[j,4]<-as.matrix(summcph$logtest)[3]
  cox_NES[j,5]<-as.matrix(summcph$sctest)[3]
  cox_NES[j,6]<-summcph$coefficients[5]
}
rownames(cox_NES)=colnames(data2)[1:(ncol(data2)-2)]
colnames(cox_NES)<-c("HR","Lower.95","Upper.95","Logtest","Logrank","p_value")

datac <- cox_NES
genes<-rownames(datac)
count<-rep("N=939", 5)
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
datac$p[datac$p_value >  0.05]<-""
datac$p.value<-paste(datac$p_value,datac$p,sep = "  ")
datac$Hazard_ratio<-paste(sprintf("%0.3f", datac$HR), "(", 
                          sprintf("%0.3f", datac$Lower.95), "-", 
                          sprintf("%0.3f", datac$Upper.95), ")", sep = "")


## The rest of the columns in the table.
tabletextc <- cbind(c("",datac$genes),
                    c("Hazard ratio(95% CI)",datac$Hazard_ratio),
                    c("P Value",datac$p.value))
library(forestplot)
## Mapping of forests

#------------------------------------------------------------
#    Figure 7A
#------------------------------------------------------------

forestplot(labeltext=tabletextc,
           graph.pos=2, #Where the Pvalue box line chart is located
           mean=c(NA,datac$HR),
           lower=c(NA,datac$Lower.95), 
           upper=c(NA,datac$Upper.95),
           #Define Title
           title="",
           ##Define the x-axis
           #xlab="    <---PCI Better---   ---Medical Therapy Better--->",
           ##Setting the line pattern according to the position of the subgroups, the width creates a "blockiness".
           
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                           "2" = gpar(lty=2),
                           #"4" = gpar(lwd=30, lineend="butt", columns=c(1:4), col="#99999922"),
                           "7"= gpar(lty=2) ),
           #The cex parameter in the fpTxtGp function sets the size of each component
           
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           ##The fpColors function sets the colors
           col=fpColors(box = "red", lines="black", zero = "gray50"),
           #Location of baseline in box plot
           zero=1,
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #Box size, line width
           lwd.ci=1, boxsize=0.3,
           #Add small vertical lines at the ends of the box-and-line diagram, height
           ci.vertices=F, ci.vertices.height = 0.4)

#-----------------------------------------------------------
#-----------------------------------------------------------

ddist <- datadist(data3)
options(datadist="ddist")
data4<-data3
f1 <- cph(Surv(OS.time, OS) ~ Age+`Risk Score`+Stage, 
          data = data4 ,surv=T,x=TRUE, y=TRUE,time.inc=365)
P1 <- psm(Surv(OS.time, OS) ~ Age+`Risk Score`+Stage, data = data3,
          x=T, y=T,dist='weibull')
data5<-data3
colnames(data5)[8]<-"RiskScore"
f2 <- cph(Surv(OS.time, OS) ~ Age+RiskScore+Stage, 
          data = data5 ,surv=T,x=TRUE, y=TRUE,time.inc=365)




#------------------------------------------------------------
#    Figure 7C
#------------------------------------------------------------
## nomogram
surv <- Survival(f1)
nom <- nomogram(f1, 
                fun=list(function(x) surv(3*365, x), 
                         function(x) surv(5*365, x),
                         function(x) surv(7*365, x), 
                         function(x) surv(10*365, x)), 
                funlabel=c("3-year survival", 
                           "5-year survival", 
                           "7-year survival", 
                           "10-year survival"),
                maxscale=100, fun.at=c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99))
plot(nom)
#-----------------------------------------------------------
#-----------------------------------------------------------

## Plotting the calibration curve
cal3 <- calibrate(P1, cmethod='KM', method='boot', u=3*365,# "u" needs to be consistent with the "time.inc" defined in the previous model
                  m=200, B=50)
cal5 <- calibrate(P1, cmethod='KM', method='boot', u=5*365,# "u" needs to be consistent with the "time.inc" defined in the previous model
                  m=200, B=50)
cal7 <- calibrate(P1, cmethod='KM', method='boot', u=7*365,# "u" needs to be consistent with the "time.inc" defined in the previous model
                  m=200, B=50)
cal10<- calibrate(P1, cmethod='KM', method="boot",u=10*365,# "u" needs to be consistent with the "time.inc" defined in the previous model
                  m=200, B=50)
#------------------------------------------------------------
#    Figure 7D
#------------------------------------------------------------

par(mar=c(5,5,2,2)) 
plot(cal3 ,lwd=2,lty=1,errbar.col='green',xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability",ylab="Actual Survival Proportion",
     subtitles = F, col='green',cex.axis=1.2,cex.lab=1.2,
     mgp=c(2,1,0))
plot(cal5 ,add = T,lwd=2,lty=1,errbar.col='yellow',xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability",ylab="Actual Survival Proportion",
     subtitles = F, col='yellow',cex.axis=1.2,cex.lab=1.2,
     mgp=c(2,1,0))
plot(cal7 ,add = T,lwd=2,lty=1,errbar.col='red',xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability",ylab="Actual Survival Proportion",
     subtitles = F, col='red',cex.axis=1.2,cex.lab=1.2,
     mgp=c(2,1,0))
plot(cal10,add = T,lwd=2,lty=1,errbar.col='blue',xlim=c(0,1),ylim=c(0,1),
     xlab="Nomogram-Predicted Probability",ylab="Actual Survival Proportion",
     subtitles = F, col='blue',cex.axis=1.2,cex.lab=1.2,
     mgp=c(2,1,0))
title(main='Training OS')
legend("topleft",
       bty='n',
       title = NULL,
       c("3-Year OS","5-Years OS","7-Years OS","10-Years OS"),
       lty = 1,
       pch = 16,
       col=c("green","yellow","red","blue"))

#-----------------------------------------------------------
#-----------------------------------------------------------

### Plotting decision curves
d_train <- dca(f2,times=c(3*365,5*365,7*365,10*365))
d_train$model<-as.vector(d_train$model)
table(d_train$model)
class(d_train$model)

d_train$model[d_train$model %in% "f2-1095"] <- "3-year"
d_train$model[d_train$model %in% "f2-1825"] <- "5-year"
d_train$model[d_train$model %in% "f2-2555"] <- "7-year"
d_train$model[d_train$model %in% "f2-3650"] <- "10-year"
d_train$model[d_train$model %in% "All-1095"] <- "3-year all"
d_train$model[d_train$model %in% "All-1825"] <- "5-year all"
d_train$model[d_train$model %in% "All-2555"] <- "7-year all"
d_train$model[d_train$model %in% "All-3650"] <- "10-year all"
d_train1<-d_train
d_train1$model<-factor(d_train1$model,levels = c("3-year","3-year all",
                                                 "5-year","5-year all",
                                                 "7-year", "7-year all",
                                                 "10-year","10-year all",
                                                 "None"))
#------------------------------------------------------------
#    Figure 7E
#------------------------------------------------------------

ggplot(d_train1,
       color = c("green","green","yellow","yellow","red","red","blue","blue","black"),
       linetype=c(1,2,1,2,1,2,1,2,1))+ 
  xlab("Threshold probability")+
  ggtitle("Test OS")+
  theme_par()+
  theme(axis.text = element_text(size=24),
        legend.title=element_blank(),
        axis.title = element_text(size=26),
        title = element_text(size=24),
        plot.title = element_text(hjust=0.5),
        legend.text = element_text(size=18),
        legend.position = c(0.8,0.7),
        legend.background = element_rect(fill = NA))

#-----------------------------------------------------------
#-----------------------------------------------------------

########################## Data processing before cox
library(survminer)
library(survival)
library(survivalROC)
library(rms)
library(tableone)
library(ggthemes)

setwd("")
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
clin<-read.table("TCGA_CDR_clinical_data.txt",header = T,sep = "\t")
os<-clin[,colnames(clin) %in% "OS" | colnames(clin) %in% "OS.time"]
rownames(os)=clin$bcr_patient_barcod
os=na.omit(os[-which(os$OS.time==0),])
expr<-training[,-c(1,2)]
name1=substr(rownames(expr),1,12)
rownames(expr)=name1
os1=os[intersect(name1,rownames(os)),]
expr=expr[intersect(name1,rownames(os)),]
dat=cbind(expr,os1)

#############################COX
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
os<-clin[,colnames(clin) %in% "OS" | 
           colnames(clin) %in% "OS.time"| 
           colnames(clin) %in% "ajcc_pathologic_tumor_stage"| 
           colnames(clin) %in% "age_at_initial_pathologic_diagnosis" ]

rownames(os)=clin$bcr_patient_barcod
os=na.omit(os[-which(os$OS.time==0),])

wgd<-training[,c(1,2)]
name1=substr(rownames(wgd),1,12)
rownames(wgd)=name1
os1=os[intersect(name1,rownames(os)),]
wgd=wgd[intersect(name1,rownames(os)),]
data=cbind(wgd,os1)
data$patient_barcode <- substr(data$array,1,12)
data1<-merge(data,risk_score,by = "patient_barcode",)
colnames(data1)<-c("patient_barcode","array","WGD",
                   "age","tumor_stage","OS",
                   "OS.time","risk_score","group")

table(data1$tumor_stage,data1$OS)
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
                   "OS","OS.time","Risk Score","Group","SubType")


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

## nomogram
f1 <- cph(Surv(OS.time, OS) ~ Age+`Risk Score`+Stage, 
          data = data3 ,surv=T,x=TRUE, y=TRUE,time.inc=365)
ddist <- datadist(data3)
surv <- Survival(f1)
nom <- nomogram(f1, 
                fun=list(function(x) surv(3*365, x), 
                         function(x) surv(5*365, x),
                         function(x) surv(7*365, x), 
                         function(x) surv(10*365, x)), 
                funlabel=c("3-year survival", 
                           "5-year survival", 
                           "7-year survival", 
                           "10-year survival"),
                maxscale=100, fun.at=c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99))

#Age
x<-nom$Age$Age
y<-nom$Age$points
a1<-lsfit(x,y)
#Risk Score
x<-nom$`Risk Score`$`Risk Score`
y<-nom$`Risk Score`$points
a2<-lsfit(x,y)
#Stage
x<-nom$Stage$Stage
y<-nom$Stage$points
a3<-lsfit(x,y)

data3$normogram <- a1$coefficients[2] *data3$Age + a1$coefficients[1]  +
  a2$coefficients[2] *data3$`Risk Score` + a2$coefficients[1]  +
  a3$coefficients[2] *data3$Stage + a3$coefficients[1]


#------------------------------------------------------------
#    Figure 7F NO.1
#------------------------------------------------------------

nobs <- NROW(data3)
surv_roc1=survivalROC(data3$OS.time, data3$OS, data3$`Risk Score`, predict.time=365*3,  method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc1$FP, surv_roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="seagreen2" , lwd=2,
     cex.axis  = 1.5,cex.lab = 1.5,cex.main = 1.5,
     #main = "Training OS",
     #main = "Test OS",
     main = "Entire OS",
     xlab="1-Specificity",ylab="Sensitivity")

par(new=TRUE)
surv_roc2=survivalROC(data3$OS.time, data3$OS, data3$Age, predict.time=365*3, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc2$FP, surv_roc2$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="deepskyblue" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc3=survivalROC(data3$OS.time, data3$OS, data3$Stage, predict.time=365*3, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc3$FP, surv_roc3$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="yellow" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")

par(new=TRUE)
surv_roc4=survivalROC(data3$OS.time, data3$OS, data3$normogram, predict.time=365*3, method = "NNE",span = 0.25*nobs^(-0.20))
plot(surv_roc4$FP, surv_roc4$TP, type="l", xlim=c(0,1), ylim=c(0,1), col="red" , lwd=2,
     cex.axis  = 1.5,cex.lab   = 1.5,
     xlab="",ylab="")

lines(x=-1:2,y=-1:2,lwd=2,col="grey60")


auc1=round(surv_roc1$AUC,2)
auc2=round(surv_roc2$AUC,2)
auc3=round(surv_roc3$AUC,2)
auc4=round(surv_roc4$AUC,2)

legend("topleft",inset=c(0.3,0.5),
       c("AUC of 3-year",
         paste("Risk Score (",auc1,")",sep = ""),
         paste("Age (",auc2,")",sep = ""),
         paste("Stage (",auc3,")",sep = ""),
         paste("Normogram (",auc4,")",sep = "")),
       col=c(NA,"seagreen2", "deepskyblue","yellow","red"),
       cex = 1.1,lty =1,lwd = 2, bty = "n")
#-----------------------------------------------------------
#-----------------------------------------------------------
