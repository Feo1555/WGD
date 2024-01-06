library(clusterProfiler)
library(stringr) 
library(AnnotationDbi)
library(org.Hs.eg.db) 
library(DOSE)
library(ggplot2) 
library(ggrepel)
setwd("")
###################################  GO
overlap<-read.csv("BRCA red&blue_module & linear_module overlap_genes_247genes.csv")
overlap<-overlap[,1]
id_list <- mapIds(org.Hs.eg.db,overlap,"ENTREZID","SYMBOL")
# Remove the result of an unsuccessful conversion, i.e., id=NA
id_list <- na.omit(id_list)
go <- enrichGO(gene = id_list, # Entrez ID list
               OrgDb = org.Hs.eg.db, # Specified species database
               keyType = "ENTREZID", # Specify the given name type
               ont = "ALL", # Optional, BP (biological process)/CC (cellular component)/MF (molecular function)/ALL (simultaneous specification)
               pAdjustMethod = "fdr", # P-value correction methods, also BH
               pvalueCutoff = 0.05,qvalueCutoff = 0.2, # p/q value threshold
               readable = T # Convert ID to symbol
)
go.res <- data.frame(go) 
go.res$Description<-gsub("[ \n ]", " ", go.res$Description)
#write.csv(go.res,"Table_GO_result.csv") # 输出GO富集分析结果

# Plot the GO enrichment analysis bar graph, the results default by qvalue ascending order, respectively, selected the top ten term to plot can be
goBP <- subset(go.res,subset = (ONTOLOGY == "BP"))
goCC <- subset(go.res,subset = (ONTOLOGY == "CC"))
goMF <- subset(go.res,subset = (ONTOLOGY == "MF"))
go.df <- rbind(goBP[1:10,],goCC[1:10,],goMF[1:10,])
#go.df<-na.omit(go.df)
# Make the order of the GO terms drawn consistent with the inputs
goBP$Description <- factor(goBP$Description,levels = rev(goBP$Description))
goCC$Description <- factor(goCC$Description,levels = rev(goCC$Description))
goMF$Description <- factor(goMF$Description,levels = rev(goMF$Description))
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
# 绘图

#------------------------------------------------------------
#    Supplementary Figure 11A
#------------------------------------------------------------

go.df$number <- factor(go.df$Description,levels = rev(go.df$Description))
ggplot(data = go.df, # Data used for mapping
       aes(x = Count,y = reorder(number,Count)))+ # Horizontal and vertical coordinates and ordering
  geom_point(aes(size = Count,color = p.adjust))+ # Bubble size and color settings
  theme_bw()+ # Remove background color
  scale_colour_gradient(low = "red",high = "purple")+ # Setting the bubble gradient color
  labs(x = "", y = "",title = "", # Setting Axis Titles and Plot Titles
       color = expression(p.adjust),size = "Count")+ # Setting the legend color and size
  guides(color = guide_colorbar(order = 1),
         size = guide_legend(order = 0))+# Setting the legend order
  facet_grid(ONTOLOGY ~ ., scale="free")+
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))

#-----------------------------------------------------------
#-----------------------------------------------------------

#######################################################################  KEGG
kegg <- enrichKEGG(gene = id_list, 
                   organism = "hsa",keyType = "kegg", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.2,
                   minGSSize = 10,maxGSSize = 500,use_internal_data = F)
# Convert IDs in the result table to symbols (not required, depending on individual needs)
kk <- setReadable(kegg,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
kegg.df <- data.frame(kk) # 结果转化成数据框
#write.csv(kegg.df,"KEGG_Table.csv")
kegg.df1$number <- factor(kegg.df1$Description,levels = rev(kegg.df1$Description))

#------------------------------------------------------------
#    Supplementary Figure 11B
#------------------------------------------------------------

ggplot(data = kegg.df1, # Data used for mapping
       aes(x = Count,y = reorder(number,Count)))+ # Horizontal and vertical coordinates and ordering
  geom_point(aes(size = Count,color = p.adjust))+ # Bubble size and color settings
  theme_bw()+ # Remove background color
  scale_colour_gradient(low = "red",high = "purple")+ # Setting the bubble gradient color
  labs(x = "", y = "",title = "", # Setting Axis Titles and Plot Titles
       color = expression(p.adjust),size = "Count")+ # Setting the legend color and size
  guides(color = guide_colorbar(order = 1),
         size = guide_legend(order = 0))+# Setting the legend order
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))

#-----------------------------------------------------------
#-----------------------------------------------------------
