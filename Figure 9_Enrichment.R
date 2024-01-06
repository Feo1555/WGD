library(ggplot2)
library(ggrepel)
library(limma)
library(GEOquery)
library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)  
library(ggupset)
library(annotate)
library(R.utils)
library(ggpubr)
library(GOplot)
library(stringr)
setwd("")
#  "BRCA_expr.csv"includes samples from the TCGA BRCA dataset with WGDs of 0 and 1, 
# and deletes genes with expression values of 0 in more than 90% of the samples.
BRCA<-read.csv("BRCA_expr.csv",row.names = 1)
clin_rs_train<-read.csv("clin_riskscore_train.csv",row.names = 1)
clin_rs_train<-clin_rs_train[,c(2,10)]
brca_rs<-merge(BRCA,clin_rs_train,by = "array")
brca_rs$group_risk

d<-brca_rs
Low<-subset(d,group_risk == "Low")
High<-subset(d,group_risk == "High")
Expression<-rbind (Low,High)
Expression<-Expression[,-c(1:3,19419)]
Expression<-t(Expression)
group_list<-c(rep("Low", nrow(Low)), rep("High", nrow(High)))
group_list<-factor(group_list,levels=c("Low" , "High"))
design<-model.matrix(~group_list)
norm<-Expression

fit<-lmFit(norm,design)
fit2<-eBayes (fit)
DEGs=topTable(fit2,coef=2,adjust.method='fdr',n=Inf,sort.by="logFC")

DEGs_new<-DEGs[rev(order(DEGs[,1])),]
DEGs_new2<-cbind(rownames(DEGs_new) ,DEGs_new$logFC, DEGs_new$adj.P.Val)
colnames(DEGs_new2)<-c("genes", "logFC", "adj")
DEGs_new2<-as.data.frame (DEGs_new2)
DEGs_new2[,2]<-as.numeric(DEGs_new2[,2])
DEGs_new2[,3]<-as.numeric(DEGs_new2[,3])
brca_DEG<-subset(DEGs_new2, abs(logFC)>0.5 & adj < 0.01)
dataset <- DEGs_new

# Setting the pvalue and logFC thresholds
cut_off_pvalue = 0.01
cut_off_logFC = 0.5
dataset$change = ifelse(dataset$P.Value < cut_off_pvalue & abs(dataset$logFC) >= cut_off_logFC, 
                        ifelse(dataset$logFC> cut_off_logFC ,'Up','Down'),
                        'Stable')
cut_off_adjpvalue = 0.01
cut_off_logFC = 0.5
table(dataset$change)

#------------------------------------------------------------
#    Figure 9A
#------------------------------------------------------------

##### volcanic composition
ggplot(dataset, 
       aes(x = logFC, 
           y = -log10(adj.P.Val), 
           colour=change)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("#4393C3", "#d2dae2","#D6604D"))+
  
  # polyline
  geom_vline(xintercept=c(-0.5,0.5),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_adjpvalue),lty=4,col="black",lwd=0.8) +
  
  # coordinate axis
  labs(x="log2(Fold Change)",
       y="-log10(FDR)")+
  theme_bw()+
  # legend 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 13),
        legend.title = element_blank()
  ) 

#-----------------------------------------------------------
#-----------------------------------------------------------

##### GSEA KEGG
#### The transformed gene is called entrez ID ####
options(stringsAsFactors = F)


need_DEG <- dataset
need_DEG <- need_DEG[,c(1,7)]   
colnames(need_DEG) <- c('log2FoldChange','change')
need_DEG$SYMBOL <- rownames(need_DEG)

##### Create geneList for gsea analysis (contains log2FoldChange and ENTREZID information sorted from largest to smallest)####
#Conversion id  
df <- bitr(rownames(need_DEG), 
           fromType = "SYMBOL",
           toType =  "ENTREZID",
           OrgDb = "org.Hs.eg.db") #Human database "org.Hs.eg.db" Mouse "org.Mm.eg.db"
need_DEG <- merge(need_DEG, df, by='SYMBOL')  #Merge comment information according to SYMBOL
geneList <- need_DEG$log2FoldChange
names(geneList) <- need_DEG$ENTREZID
geneList <- sort(geneList, decreasing = T)   #sort
R.utils::setOption("clusterProfiler.download.method",'auto')
KEGG_kk_entrez <- gseKEGG(geneList = geneList,
                          organism = "hsa", #human "hsa" rat "mmu"
                          pvalueCutoff = 0.05)  #Actual for padj threshold adjustable 
KEGG_kk <- DOSE::setReadable(KEGG_kk_entrez,
                             OrgDb="org.Hs.eg.db",
                             keyType='ENTREZID')# Conversion id

p<-KEGG_kk_entrez@result$ID[KEGG_kk_entrez@result$NES > 0]
px<-KEGG_kk_entrez@result$ID[order(KEGG_kk_entrez@result$pvalue,decreasing=FALSE)[1:10]]
KEGG_kk_entrez@result$Description<-gsub("\\(.*?\\)","",gsub(" - Homo sapiens ","",KEGG_kk_entrez@result$Description))
#------------------------------------------------------------
#    Figure 9B NO.1 
#------------------------------------------------------------

gseaplot2(KEGG_kk_entrez, 
          p,  
          rel_heights = c(1.8, 0.5, 0.7),#Relative height of the sub-chart
          base_size = 10,
          pvalue_table = F)

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Figure 9B NO.2 
#------------------------------------------------------------

gseaplot2(KEGG_kk_entrez, 
          px,  
          rel_heights = c(1.8, 0.5, 0.7),#Relative height of the sub-chart
          base_size = 10,
          pvalue_table = F)

#-----------------------------------------------------------
#-----------------------------------------------------------


#####  Filter settings #######
log2FC_cutoff = 0.5
pvalue_cutoff = 0.01
padj_cutoff = 0.01

####Obtaining up- and down-regulated genes for DEG results ####

DEG<-dataset
need_DEG <- DEG[,c(1,4,5)]  
head(need_DEG)
colnames(need_DEG) <- c('log2FoldChange','pvalue','padj')
rownames(dataset)[dataset$change %in% "Up"]
gene_up=rownames(dataset)[dataset$change %in% "Up"]
gene_down=rownames(dataset)[dataset$change %in% "Down"]

#### The transformed gene is called entrez ID ####
#org.Hs.eg.db\org.Mm.eg.db contains data from all the major databases, such as entrez ID and ensembl,
#keytypes(org.Hs.eg.db) #View all supported and transformable types Commonly "ENTREZID", "ENSEMBL", "SYMBOL"
gene_up_entrez <- as.character(na.omit(bitr(gene_up,
                                            fromType="SYMBOL", #input format
                                            toType="ENTREZID", # Conversion to ENTERZID format
                                            OrgDb="org.Hs.eg.db")[,2])) #"org.Hs.eg.db" "org.Mm.eg.db"
gene_down_entrez <- as.character(na.omit(bitr(gene_down, 
                                              fromType="SYMBOL", #input format
                                              toType="ENTREZID", # Conversion to ENTERZID format
                                              OrgDb="org.Hs.eg.db")[,2])) #"org.Hs.eg.db" "org.Mm.eg.db"
gene_diff_entrez <- unique(c(gene_up_entrez ,gene_down_entrez ))


#############  KEGG
R.utils::setOption("clusterProfiler.download.method",'auto')
kegg_up_results <- enrichKEGG(gene  = gene_up_entrez,
                              organism  = "hsa",
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)
kk_up <- DOSE::setReadable(kegg_up_results,
                           OrgDb="org.Hs.eg.db",
                           keyType='ENTREZID')#ENTREZID to gene Symbol
#write.csv(kk_up@result,'KEGG_gene_up.csv')


kegg_down_results <- enrichKEGG(gene  = gene_up_entrez,
                                organism  = "hsa",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2)
kk_down <- DOSE::setReadable(kegg_down_results,
                             OrgDb="org.Hs.eg.db",
                             keyType='ENTREZID')#ENTREZID to gene Symbol
#write.csv(kk_down@result,'KEGG_gene_down.csv')
info<-dataset
info$SYMBOL<-rownames(info)
genedata<-data.frame(ID=info$SYMBOL,logFC=info$logFC)

############################   KEGG
kegg_up <- kk_up@result
kegg_up<-subset(kegg_up,p.adjust <0.05)
keggupplotIn<-kegg_up[,c("ID", "Description", "p.adjust", "geneID")] #Extract the first 10 rows of kegg-enriched BPs, and extract the four columns of ID, Description, p.adjust, and GeneID.
keggupplotIn$geneID <-str_replace_all(keggupplotIn$geneID,'/',',')
names(keggupplotIn)<-c('ID','Term','adj_pval','Genes')#Modify the column names, this format is needed for the later string diagrams.
keggupplotIn$category="kegg"
circ_keggup<-GOplot::circle_dat(keggupplotIn,genedata) #goplot import data format organization

kegg_down <- kk_down@result
kegg_down<-subset(kegg_down,p.adjust <0.05)
keggdownplotIn<-kegg_down[,c("ID", "Description", "p.adjust", "geneID")] 
keggdownplotIn$geneID <-str_replace_all(keggdownplotIn$geneID,'/',',')
names(keggdownplotIn)<-c('ID','Term','adj_pval','Genes')
keggdownplotIn$category="kegg"
circ_keggdown<-GOplot::circle_dat(keggdownplotIn,genedata) 


# Plotting
# reorder so that the vertical axis is sorted by Count

#------------------------------------------------------------
#    Figure 9A NO.1
#------------------------------------------------------------

p1<-GOCircle(circ_keggup,
             title = "Upregulated in High Risk",
             rad1 = 2, rad2 = 3, # I.D., O.D. Setting
             lfc.col = 'firebrick3',# Upward Color Settings
             label.size = 5,
             label.fontface='bold',# Font Size Formatting
             table.legend = F 
             # The right side of the table is set, TRUE prevents the theme from being set.
             # If you unset the table to add the theme, it will release the gridlines and axes.
) +
  guides(color = guide_colorbar(order = 1),
         size = guide_legend(order = 0))+# Setting the legend order
  # Remove large background color
  theme_bw() +
  theme(
    # Remove axes text, tick marks, titles
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 20),
    # Remove axes borders, gridlines in graphs
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  # Font-related settings
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom" )

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Figure 9C NO.2
#------------------------------------------------------------

p2<-GOCircle(circ_keggdown,
             title = "Downregulated in High Risk",
             rad1 = 2, rad2 = 3, 
             lfc.col = 'royalblue3',
             label.size = 5,
             label.fontface='bold',
             table.legend = F 
) +
  theme_bw() +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    title = element_text(size = 20),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.position = "bottom" )

#-----------------------------------------------------------
#-----------------------------------------------------------

############  GO
go_up <- enrichGO(gene = gene_up_entrez,
                  OrgDb = "org.Hs.eg.db",
                  ont   = "ALL"  ,     #One of "BP", "MF"  "CC"  "ALL"
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
#write.csv(go_up@result, 'GO_gene_up.csv')

go_down <- enrichGO(gene = gene_down_entrez,
                    OrgDb = "org.Hs.eg.db",
                    ont   = "ALL"  ,     #One of "BP", "MF"  "CC"  "ALL"
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    readable      = TRUE)
#write.csv(go_down@result, 'GO_gene_down.csv')


# Plot the GO enrichment analysis bar graph, the results default by qvalue ascending order, respectively, selected the top ten term to plot can be
go_up<-go_up@result
goBP <- subset(go_up,subset = (ONTOLOGY == "BP"))

# Make the order of the GO terms drawn consistent with the inputs
goBP$Description <- factor(goBP$Description,levels = rev(goBP$Description))

# plot
goBP$number <- factor(goBP$Description,levels = rev(goBP$Description))
up<-goBP

#------------------------------------------------------------
#    Figure 9D NO.1
#------------------------------------------------------------

p1<-ggplot(data = up, 
           aes(x = Count,y = reorder(number,Count)))+ 
  geom_point(aes(size = Count,color = p.adjust))+ 
  theme_bw()+ # 去除背景色
  scale_colour_gradient(low = "red",high = "purple")+ 
  labs(x = "", y = "",title = "Upregulated in High Risk", 
       color = expression(p.adjust),size = "Count")+
  guides(color = guide_colorbar(order = 1),
         size = guide_legend(order = 0))+
  facet_grid(ONTOLOGY ~ ., scale="free")+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10.5))

#-----------------------------------------------------------
#-----------------------------------------------------------

go_down<-go_down@result
goBP <- subset(go_down,subset = (ONTOLOGY == "BP"))
goCC <- subset(go_down,subset = (ONTOLOGY == "CC"))
goMF <- subset(go_down,subset = (ONTOLOGY == "MF"))
go.df <- rbind(goBP[1:10,],goCC[1:10,],goMF[1:10,])
go.df<-na.omit(go.df)
goBP$Description <- factor(goBP$Description,levels = rev(goBP$Description))
goCC$Description <- factor(goCC$Description,levels = rev(goCC$Description))
goMF$Description <- factor(goMF$Description,levels = rev(goMF$Description))
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
# plot
go.df$number <- factor(go.df$Description,levels = rev(go.df$Description))

#------------------------------------------------------------
#    Figure 9D NO.2
#------------------------------------------------------------
p2<-ggplot(data = go.df, 
           aes(x = Count,y = reorder(number,Count)))+ 
  geom_point(aes(size = Count,color = p.adjust))+ 
  theme_bw()+ 
  scale_colour_gradient(low = "red",high = "purple")+ 
  labs(x = "", y = "",title = "Downregulated in Low Risk", 
       color = expression(p.adjust),size = "Count")+ 
  guides(color = guide_colorbar(order = 1),
         size = guide_legend(order = 0))+
  facet_grid(ONTOLOGY ~ ., scale="free")+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10.5))

#-----------------------------------------------------------
#-----------------------------------------------------------