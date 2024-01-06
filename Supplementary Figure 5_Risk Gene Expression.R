library(ggplot2)
#Defining Graphical Functions 
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL){
                             # By @YAK: https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1,'group']
                             newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           }
)
create_quantile_segment_frame <- function (data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density)/sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(x = ggplot2:::interleave(violin.xs, violin.xmaxvs), 
               y = rep(ys, each = 2), group = rep(ys, each = 2)) 
  } else {
    data.frame(x = ggplot2:::interleave(violin.xminvs, violin.xs), 
               y = rep(ys, each = 2), group = rep(ys, each = 2)) 
  }
}
geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

###

library(reshape2)    
library(ggpubr)
#"22gene_exp_train.RData" is the expression of 22gene in different risk grouping samples of the training cohort
load("22gene_exp_train.RData")

#------------------------------------------------------------
#    Supplementary Figure 5 NO.1
#------------------------------------------------------------

p1<-ggplot(data1, aes(x =Genes, y = `Gene expression`, fill = group)) +
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75),lwd = 0.7) + 
  xlab("")+
  stat_compare_means(aes(group=group),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", " ")),
                     label = "p.signif")+
  scale_fill_manual(values = c("#FF8C69","#87CEFF"))+
  scale_color_manual(values = "black")+
  ggtitle("Training Cohort")+
  theme_classic()+
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        #panel.background = element_rect(fill = "#F0FFF0",colour = "#F0FFF0"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        #text = element_text(size = 10),
        legend.position = 'top',
        panel.border = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size=20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle=90, hjust=1, vjust=1)
        #legend.key.size = unit(0.3, "inches")
  )

#-----------------------------------------------------------
#-----------------------------------------------------------

#"22gene_exp_test.RData" is the expression of 22gene in different risk grouping samples of the test cohort 
load("22gene_exp_test.RData")

#------------------------------------------------------------
#    Supplementary Figure 5 NO.2
#------------------------------------------------------------

p2<-ggplot(data1, aes(x =Genes, y = `Gene expression`, fill = group)) +
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75),lwd = 0.7) + 
  scale_fill_manual(values = c("#FF8C69","#87CEFF"))+
  scale_color_manual(values = c("#FF8C69","#87CEFF"))+
  stat_compare_means(aes(group=group),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", " ")),
                     label = "p.signif")+
  xlab("")+
  ggtitle("Test Cohort")+
  theme_classic()+
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        #panel.background = element_rect(fill = "#F0FFF0",colour = "#F0FFF0"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        #text = element_text(size = 10),
        legend.position = 'top',
        panel.border = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size=20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle=90, hjust=1, vjust=1)
        #legend.key.size = unit(0.3, "inches")
  )

#-----------------------------------------------------------
#-----------------------------------------------------------

#"22gene_exp_entire.RData" is the expression of 22gene in different risk grouping samples of the entire cohort 
load("22gene_exp_entire.RData")

#------------------------------------------------------------
#    Supplementary Figure 5 NO.3
#------------------------------------------------------------

p3<-ggplot(data1, aes(x =Genes, y = `Gene expression`, fill = group)) +
  geom_split_violin(draw_quantiles = c(0.25, 0.5, 0.75),lwd = 0.7) + 
  xlab("")+
  stat_compare_means(aes(group=group),
                     method="wilcox.test",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                      symbols = c("***", "**", "*", " ")),
                     label = "p.signif")+
  scale_fill_manual(values = c("#FF8C69","#87CEFF"))+
  scale_color_manual(values = "black")+
  ggtitle("Entire Cohort")+
  theme_classic()+
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        #panel.grid.major=element_line(colour=NA),
        #panel.background = element_rect(fill = "#F0FFF0",colour = "#F0FFF0"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        #text = element_text(size = 10),
        legend.position = 'top',
        panel.border = element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size=20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle=90, hjust=1, vjust=1)
        #legend.key.size = unit(0.3, "inches")
  )

#-----------------------------------------------------------
#-----------------------------------------------------------
