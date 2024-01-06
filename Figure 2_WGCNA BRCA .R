library(WGCNA)
library(AnnoProbe)
library(GEOquery)
library(stringr)
library(readr)
library(pheatmap)
library(ggplot2)
library(openxlsx)
library(pheatmap)
library(limma)
library(ggpubr)
library(ggthemes)

# upir_num_all_log1.Rdata is RNA expression data for BRCA from TCGA, 
# its rows are samples columns are genes
load("upir_num_all_log1.Rdata")
####################################################################ploid,purity,WGD
ploidy_f<-read.table("sample-purity-ploidy-WGD-cancer.type.txt")
ploidy<-cbind(ploidy_f$array,
              ploidy_f$ploidy,
              ploidy_f$purity,
              ploidy_f$Genome.doublings)
colnames(ploidy)<-c("array","ploidy","purity","WGD")
table(upir_num_all_log1$cancer.type)
a<-unique(upir_num_all_log1$cancer.type)
disease<-a[c(3,4,7,15)]

########################################################## variance
data0<-subset(upir_num_all_log1,!WGD %in% 2)
#去掉patient，wgd，cancer.type
data<-data0[,-c(1,ncol(data0)-1,ncol(data0))]

variance<-as.data.frame(apply(data, 2, var))
variance[,2]<-rownames(variance)
variance[,3]<-variance[,1]
variance<-variance[,-1]
colnames(variance)<-c("gene","var")
var<-variance[order(variance[,2],decreasing=T),]
c<-nrow(var)*0.25
var_gene<-var[1:c,1]

#####################################################  WGCNA

expr<-subset(upir_num_all_log1,cancer.type %in% disease[i]&WGD != "2")
expr_var<-expr[,colnames(expr) %in% var_gene]
datExpr0<-expr_var
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK
#logistic
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


meanFPKM=0.5  
n=nrow(datExpr0)
datExpr0[n+1,]=apply(datExpr0[c(1:nrow(datExpr0)),],2,mean)

# for meanFpkm in row n+1 and it must be above what you set--select meanFpkm>opt$meanFpkm(by rp)
datExpr0=datExpr0[1:n,datExpr0[n+1,] > meanFPKM]

filtered_fpkm=t(datExpr0)
filtered_fpkm=data.frame(rownames(filtered_fpkm),filtered_fpkm)
names(filtered_fpkm)[1]="sample"
head(filtered_fpkm)
#write.table(filtered_fpkm, file="rt_filter.xls",row.names=F, col.names=T,quote=FALSE,sep="\t")


sampleTree = hclust(dist(datExpr0), method = "average")
samwgd = data0$WGD
table(samwgd)
colors = numbers2colors(samwgd)
table(colors)

enableWGCNAThreads() 
powers = c(1:20) 
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5,RsquaredCut = 0.9)
cex1<-0.9

#------------------------------------------------------------
#    Figure 2A
#------------------------------------------------------------

par(mfrow = c(1,2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.9,col="red") 

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Figure 2B
#------------------------------------------------------------

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#-----------------------------------------------------------
#-----------------------------------------------------------

softPower =sft$powerEstimate
adjacency = adjacency(datExpr0, power = softPower)
#dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = softPower)

TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

geneTree = hclust(as.dist(dissTOM), method = "average");



minModuleSize = 30 
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average")

moduleColors = MEList$validColors

nSamples=nrow(datExpr0)
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
names(datExpr0)
probes = names(datExpr0)


geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
#write.table(geneInfo, file = "MM.xls",sep="\t",row.names=F)
color<-unique(geneInfo$moduleColor)
#################################################################
genetocolor<-geneInfo[,1:2]
nGenes = ncol(datExpr0)
plotTOM = dissTOM^7
diag(plotTOM) = NA
par(cex = 0.9)

#------------------------------------------------------------
#    Figure 2C
#------------------------------------------------------------

plotEigengeneNetworks(MEs, 
                      "Eigengene dendrogram & Eigengene adjacency heatmap", 
                      marDendro = c(2,3,2,1), 
                      marHeatmap = c(2,3,2,1), 
                      cex.lab = 1, 
                      xLabelsAngle= 90)

#-----------------------------------------------------------
#-----------------------------------------------------------

cytoDir="CytoscapeInput"
#dir.create(cytoDir)
for (mod in 1:nrow(table(moduleColors))){
  modules = names(table(moduleColors))[mod]
  probes = names(datExpr0)
  inModule = (moduleColors == modules)
  modProbes = probes[inModule]
  modGenes = modProbes
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  edges_File = paste("CytoscapeInput-edges-", modules , ".txt", sep="")
  nodes_File = paste("CytoscapeInput-nodes-", modules, ".txt", sep="")
  outEdge=paste(cytoDir,edges_File,sep="\\")
  outNode=paste(cytoDir,nodes_File,sep="\\")
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = outEdge,
                                 nodeFile = outNode,
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modProbes,
                                 altNodeNames = modGenes,
                                 nodeAttr = moduleColors[inModule])
}



nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

test<-matrix(0,nrow(datExpr0),2)
rownames(test)<-rownames(datExpr0)
test[,1]<-substr(rownames(test),1,15)
test[,2]<-rownames(test)
colnames(test)<-c("array","sample")

feature0<-merge(test,ploidy,by = "array",all.x = TRUE)
rownames(feature0)<-feature0$sample


############################################################################### CIN

d<-"brca"
CIN_0<-read.xlsx("Supplementary File-data2.xlsx",sheet = 1,rowNames = T)
CIN_1<-read.xlsx("Supplementary File-data2.xlsx",sheet = 2,rowNames = T)


colnames(CIN_0)<-CIN_0[1,]
colnames(CIN_1)<-CIN_1[1,]
CIN_0<-CIN_0[-1,]
CIN_1<-CIN_1[-1,]

colnames(CIN_0)[3]<-"CIN"
colnames(CIN_1)[3]<-"CIN"

cin<-rbind(CIN_0,CIN_1)
cin[,4]<-rownames(cin)
cin<-cin[,c(3,4)]
colnames(cin)[2]<-c("array")
feature_cin<-merge(feature0,cin,by = "array",all.x = TRUE)
feature_cin<-unique(feature_cin)

rownames(feature_cin)<-feature_cin$sample
feature3<-feature_cin[,-c(1:2)]

feature3<-na.omit(feature3)
traitData<-feature3
corType = "pearson"
MEs_col<-MEs
robustY = ifelse(corType=="pearson",T,F)

traitData1<-traitData[,c(1,2,4)]
MEs_col=MEs_col[rownames(traitData1),]

if (corType=="pearson") {
  modTraitCor1 = cor(MEs_col, traitData1, use = "p")
  modTraitP1 = corPvalueStudent(modTraitCor1, nSamples)
} else {
  modTraitCorP1 = bicorAndPvalue(MEs_col, traitData1, robustY=robustY)
  modTraitCor1 = modTraitCorP1$bicor
  modTraitP1   = modTraitCorP1$p
}
#####################################
MEs_col<-MEs
robustY = ifelse(corType=="spearman",T,F)
traitData2<-as.data.frame(traitData[,3])
rownames(traitData2)<-rownames(traitData)
colnames(traitData2)<-colnames(traitData)[3]
MEs_col=MEs_col[rownames(traitData2),]
class(traitData2[1,1])

if (corType=="spearman") {
  modTraitCor2 = cor(MEs_col, traitData2, use = "spearman")
  modTraitP2 = corPvalueStudent(modTraitCor2, nSamples)
} else {
  modTraitCorP2 = bicorAndPvalue(MEs_col, traitData2, robustY=robustY)
  modTraitCor2 = modTraitCorP2$bicor
  modTraitP2   = modTraitCorP2$p
}
modTraitCor<-cbind(modTraitCor1,modTraitCor2)
modTraitP<-cbind(modTraitP1,modTraitP2)
# "signif" indicates the number of decimal places to be retained
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)

#------------------------------------------------------------
#    Figure 2D
#------------------------------------------------------------

par(mar = c(3, 8, 2, 2))
labeledHeatmap(Matrix = modTraitCor, 
               xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 1, 
               ySymbols = colnames(MEs_col), 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.7, 
               zlim = c(-1,1),
               main = "Module-trait relationships")

#-----------------------------------------------------------
#-----------------------------------------------------------

sample2WGD<-read.table("D:/WGD/result/ploidy_f.txt",header = T)
colnames(sample2WGD)[6]<-"WGD"
WGD0<-subset(sample2WGD,cancer.type %in% disease[i])
WGD<-WGD0[,c(2,6)]
MEs[,ncol(MEs)+1]<-substr(rownames(MEs),1,15)
colnames(MEs)[ncol(MEs)]<-"array"

MEs[,ncol(MEs)+1]<-rownames(MEs)
colnames(MEs)[ncol(MEs)]<-"sample"

test<-merge(MEs,WGD,by = "array")
rownames(test)<-MEs$sample
test<-test[,-c(1,ncol(test)-1)]

W0_t<-subset(test,WGD == 0)
W1_t<-subset(test,WGD == 1)
test1<-rbind(W0_t,W1_t)
test2<-test1[,-which(colnames(test1) == "WGD")]
test2<-as.data.frame(apply(test2,1,scale))
rownames(test2)<-colnames(test1)[-ncol(test1)]
WGD<-c(rep("0", nrow(W0_t)), rep("1", nrow(W1_t)))
WGD<-as.data.frame(WGD)
rownames(WGD) = colnames(test2)

#############################################################################

a<-which(modTraitCor==max(modTraitCor),arr.ind=T)
modul<-substr(colnames(MEs),3,nchar(colnames(MEs)))
module = modul[c(a[1,1])]

# WGD
WGD = as.data.frame(traitData$WGD) 
names(WGD) = "WGD"

modNames = substring(names(MEs), 3)

# MM_P value
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# GS
geneTraitSignificance = as.data.frame(cor(datExpr0[rownames(datExpr0) %in% rownames(traitData1),], WGD, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                          nSamples))
names(geneTraitSignificance) = paste("GS.", names(WGD), sep="")
names(GSPvalue) = paste("p.GS.", names(WGD), sep="")


column = match(module, modNames)
moduleGenes = moduleColors==module

a<-as.data.frame(cbind(abs(geneModuleMembership[moduleGenes,column]),
                       abs(geneTraitSignificance[moduleGenes, 1])))
colnames(a)<-c("MM","GS")

#------------------------------------------------------------
#    Figure 2F
#------------------------------------------------------------

ggscatter(a, 
          x = "MM", 
          y = "GS",
          color = "red",
          size = 2,)+  
  theme_par()+ 
  geom_segment(x = -Inf, y = 0.2,aes(xend=0.5, yend=0.2),
               lty = 2,
               lwd = 1,
               color = "#EEC591") +
  geom_segment(x = 0.5, y = -Inf,aes(xend=0.5, yend=0.2),
               lty = 2,
               lwd = 1,
               color = "#EEC591")+
  stat_cor(data=a, method = "pearson",size = 5)+#R = 0.78, p < 2.6e−67 
  geom_rect(aes(xmin = 0.5,
                xmax = max(MM)+0.02,
                ymin = 0.2, 
                ymax = max(GS)+0.02),
            fill = "#EEC591", alpha = 0.01)+
  labs(x = ("Module Membership in red module"),
       y = ("Gene Significance for Cluster"))+
  ggtitle("Module Membership vs Gene Significance")+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15),
        legend.position = "none")

#-----------------------------------------------------------
#-----------------------------------------------------------

###################################################################red module genes
column = match(module, modNames)
moduleGenes = moduleColors==module
gs_mm<-matrix(0,length(abs(geneModuleMembership[moduleGenes,column])),2)
rownames(gs_mm)<-rownames(geneModuleMembership)[moduleGenes]
colnames(gs_mm)<-c("gs","mm")
gs_mm[,1]<-abs(geneTraitSignificance[moduleGenes, 1])
gs_mm[,2]<-abs(geneModuleMembership[moduleGenes,column])
write.csv(gs_mm,"gs_mm_red_WGD.csv")
genes<-rownames(gs_mm)[gs_mm[,1]>0.2 & gs_mm[,2]>0.5]
write.table(genes,"red_genes.txt",col.names = F,row.names = F)
gs_mm_red<-gs_mm
#####################################################################
b<-which(modTraitCor==min(modTraitCor),arr.ind=T)
modul<-substr(colnames(MEs),3,nchar(colnames(MEs)))
module = modul[c(b[1,1])]

# purity
purity = as.data.frame(traitData$purity) 
names(purity) = "purity"

modNames = substring(names(MEs), 3)

# MM_P value
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

# GS
geneTraitSignificance = as.data.frame(cor(datExpr0[rownames(datExpr0) %in% rownames(traitData1),], purity, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                          nSamples))
names(geneTraitSignificance) = paste("GS.", names(purity), sep="")
names(GSPvalue) = paste("p.GS.", names(purity), sep="")

column = match(module, modNames)
moduleGenes = moduleColors==module

a<-as.data.frame(cbind(abs(geneModuleMembership[moduleGenes,column]),
                       abs(geneTraitSignificance[moduleGenes, 1])))
colnames(a)<-c("MM","GS")

#------------------------------------------------------------
#    Figure 2E
#------------------------------------------------------------

ggscatter(a, 
          x = "MM", 
          y = "GS",
          color = "blue",
          size = 2,)+  
  theme_par()+ 
  geom_segment(x = -Inf, y = 0.5,aes(xend=0.5, yend=0.5),
               lty = 2,
               lwd = 1,
               color = "#EEC591") +
  geom_segment(x = 0.5, y = -Inf,aes(xend=0.5, yend=0.5),
               lty = 2,
               lwd = 1,
               color = "#EEC591")+
  stat_cor(data=a, method = "pearson",size = 5)+#R = 0.95, p < 1e−200 
  geom_rect(aes(xmin = 0.5,
                xmax = max(MM)+0.02,
                ymin = 0.5, 
                ymax = max(GS)+0.02),
            fill = "#EEC591", alpha = 0.01)+
  labs(x = ("Module Membership in blue module"),
       y = ("Gene Significance for Cluster"))+
  ggtitle("Module Membership vs Gene Significance")+
  theme(axis.line =  element_line(colour = "black"),
        strip.background = element_rect(color = "white"),
        #panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 15),
        legend.position = "none")


#-----------------------------------------------------------
#-----------------------------------------------------------



###################################################################blue module genes
column = match(module, modNames)
moduleGenes = moduleColors==module
gs_mm<-matrix(0,length(abs(geneModuleMembership[moduleGenes,column])),2)
rownames(gs_mm)<-rownames(geneModuleMembership)[moduleGenes]
colnames(gs_mm)<-c("gs","mm")
gs_mm[,1]<-abs(geneTraitSignificance[moduleGenes, 1])
gs_mm[,2]<-abs(geneModuleMembership[moduleGenes,column])
write.csv(gs_mm,"gs_mm_blue_purity.csv")
genes<-rownames(gs_mm)[gs_mm[,1]>0.35 & gs_mm[,2]>0.5]
#write.table(genes,"blue_genes.txt",col.names = F,row.names = F)
gs_mm_blue<-gs_mm



