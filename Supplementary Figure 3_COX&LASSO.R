########################## Pre-cox data processing
library(survminer)
library(survival)

# "genes" is the signature of 247 genes previously obtained
genes1<-read.csv("247genes_cnv_mutation.csv",header = T,row.names = 1)
genes<-sort(genes1[,1])
genes<-gsub("-", "\\.", genes)

# "brca_247genes_expr" is the expression data of 247 genes in BRCA
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
test<-rbind(expr0[-b,],expr1[-a,])

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

cox_NES<-matrix(nrow=ncol(dat)-2,ncol=6)
for (j in 1:(ncol(dat)-2)) {
  Bcox<-coxph(Surv(OS.time, OS)~as.numeric(dat[,j]) < as.numeric(median(dat[,j])),data=dat)
  summcph<-summary(Bcox)
  cox_NES[j,1]<-summcph$conf.int[1]
  cox_NES[j,2]<-summcph$conf.int[3]
  cox_NES[j,3]<-summcph$conf.int[4]
  cox_NES[j,4]<-as.matrix(summcph$logtest)[3]
  cox_NES[j,5]<-as.matrix(summcph$sctest)[3]
  cox_NES[j,6]<-summcph$coefficients[5]
}
rownames(cox_NES)=colnames(dat)[1:(ncol(dat)-2)]
colnames(cox_NES)<-c("HR","Lower.95","Upper.95","Logtest","Logrank","p_value")
#write.csv(cox_NES,file = "247-122gene-cox.csv")
sig_cox=cox_NES[cox_NES[,6]<0.05,]


library(glmnet)
#LASSO
opar=par(no.readonly = T)
par(family="serif",ps="15")

set.seed(100)
sig_train_data=as.matrix(dat[,rownames(sig_cox)])
TS=cbind(time=dat$OS.time,status=dat$OS)
rownames(TS)=rownames(dat)
train_glm=glmnet(sig_train_data,TS,family = "cox", alpha=1)
train_glm


#------------------------------------------------------------
#    Supplementary Figure 2A
#------------------------------------------------------------
par(family="serif",ps="15")
plot(train_glm, xvar="lambda", label=F)
#-----------------------------------------------------------
#-----------------------------------------------------------



train_cvglm=cv.glmnet(sig_train_data,TS,family = "cox", alpha=1)
                      
coe.min=coef(train_cvglm,s=train_cvglm$lambda.min)
act_index <- which(coe.min != 0)
act_coe <- coe.min[act_index]
sig_glmname=row.names(coe.min)[act_index]#Get the reduced variable
sig_glmname
beta <- as.matrix(coe.min[act_index,])
colnames(beta)<-"beta"

#write(beta,file = "22beta.txt")
#------------------------------------------------------------
#    Supplementary Figure 2B
#------------------------------------------------------------
par(family="serif",ps="15")
plot(train_cvglm,xlab="Log Lambda")
#-----------------------------------------------------------
#-----------------------------------------------------------


