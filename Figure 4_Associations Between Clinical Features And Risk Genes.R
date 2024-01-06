library(tidyverse)
library(ggsignif)
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(ggalluvial)
library(dplyr)
par(mar=c(0,3,0,3))

############################   Box

#############################################    Training
#"clinic_train.csv" includes the clinical characteristics of the training group sample
data<-read.csv("clinic_train.csv",row.names = 1, check.names = F)
data$Stage<-substr(data$Stage, 6, nchar(data$Stage))

#------------------------------------------------------------
#     Figure 4A NO.1(left to right)
#------------------------------------------------------------

p1<-ggplot(data, aes(x=Stage,
                     y=`Risk Score`,
                     color=Stage),
           palette = c("#FF8C69","#87CEFF")) + 
  geom_boxplot(width=0.5,
               lwd=0.1,
               position=position_dodge(width=0.5))+
  scale_fill_nejm()+
  ggtitle("Training",)+
  geom_jitter(width = 0.2)+
  labs(x="",
       y="Risk Score")+
  geom_signif(comparisons=list(#c("StageI","StageII"),# Add significance tags
    #c("StageI","StageIII"),
    c("StageI","StageIV"),
    #c("StageII","StageIII"),
    c("StageII","StageIV"),
    c("StageIII","StageIV")), # Choose which 2 groups you want to add labels to
    step_increase = 0.1,
    test="wilcox.test", # "t-test, comparing two groups (parametric)" = "t.test", "Wilcoxon signed rank test, comparing two groups (nonparametric)" = "wilcox.test"
    map_signif_level=T  # Label style F is a number, T is a * sign
  )+
  theme_par()+
  stat_compare_means(method = "kruskal.test",
                     size = 7,
                     label = "p",
                     # label.x = 1,
                     label.y = 1.7
  )+ #method denotes the test model used, kruskal.test is a nonparametric test and t.test is a parametric test
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        panel.background = element_rect(fill = "#F5FFFA",colour = "#F5FFFA"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size=25),
        axis.title = element_text(size=25))


#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#     Figure 4A NO.4(left to right)
#------------------------------------------------------------

p2<-ggplot(data, aes(x=SubType,
                     y=`Risk Score`,
                     color=SubType),
           palette = c("#FF8C69","#87CEFF")) + 
  geom_boxplot(width=0.5,
               lwd=0.1,
               position=position_dodge(width=0.5))+
  scale_fill_nejm()+
  ggtitle("Training",)+
  geom_jitter(width = 0.2)+
  labs(x="",
       y="Risk Score")+
  geom_signif(comparisons=list(#c("Basal","Her2"),# Add significance tags
    c("Basal","LumA"),
    c("Basal","LumB"),
    c("Basal","Normal"),
    #c("Her2","LumA"),
    c("Her2","LumB"),
    c("Her2","Normal"),
    c("LumA","LumB"),
    c("LumA","Normal"),
    c("LumB","Normal")), # Choose which 2 groups you want to add labels to
    step_increase = 0.1,
    test="wilcox.test", # "t-test, comparing two groups (parametric)" = "t.test", "Wilcoxon signed rank test, comparing two groups (nonparametric)" = "wilcox.test"
    map_signif_level=T  # Label style F is a number, T is a * sign
  )+
  theme_par()+
  stat_compare_means(method = "kruskal.test",
                     size = 7,
                     label = "p",
                     # label.x = 1,
                     label.y = 3.4
  )+ #method denotes the test model used, kruskal.test is a nonparametric test and t.test is a parametric test
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        panel.background = element_rect(fill = "#F5FFFA",colour = "#F5FFFA"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle=25, hjust=1, vjust=1))

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#     Figure 4A NO.7(left to right)
#------------------------------------------------------------

p3<-ggplot(data, aes(x=Age,
                     y=`Risk Score`,
                     color=Age),
           palette = c("#FF8C69","#87CEFF")) + 
  geom_boxplot(width=0.5,
               lwd=0.1,
               position=position_dodge(width=0.5))+
  scale_fill_nejm()+
  ggtitle("Training",)+
  geom_jitter(width = 0.2)+
  labs(x="",
       y="Risk Score")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 7,
                     label.y = .8)+ 
  theme_par()+
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        panel.background = element_rect(fill = "#F5FFFA",colour = "#F5FFFA"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25))

#############################################    Test
#"clinic_test.csv" includes the clinical characteristics of the test group sample
data<-read.csv("clinic_test.csv",row.names = 1,check.names = F)
data$Stage<-substr(data$Stage, 6, nchar(data$Stage))

#------------------------------------------------------------
#     Figure 4A NO.2(left to right)
#------------------------------------------------------------

p4<-ggplot(data, aes(x=Stage,
                     y=`Risk Score`,
                     color=Stage),
           palette = c("#FF8C69","#87CEFF")) + 
  geom_boxplot(width=0.5,
               lwd=0.1,
               position=position_dodge(width=0.5))+
  scale_fill_nejm()+
  ggtitle("Test",)+
  geom_jitter(width = 0.2)+
  labs(x="",
       y="Risk Score")+
  # geom_signif(comparisons=list(
  #   c("StageI","StageII"),# Add significance tags
  #   c("StageI","StageIII"),
  #   c("StageI","StageIV"),
  #   c("StageII","StageIII"),
  #   c("StageII","StageIV"),
  #   c("StageIII","StageIV")), # Choose which 2 groups you want to add labels to
  #   step_increase = 0.1,
  #   test="wilcox.test", # "t-test, comparing two groups (parametric)" = "t.test", "Wilcoxon signed rank test, comparing two groups (nonparametric)" = "wilcox.test"
  #   map_signif_level=T  # Label style F is a number, T is a * sign
  # )+
theme_par()+
  stat_compare_means(method = "kruskal.test",
                     size = 7,
                     label = "p",
                     # label.x = 1,
                     label.y = 0.8
  )+ #method denotes the test model used, kruskal.test is a nonparametric test and t.test is a parametric test
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        panel.background = element_rect(fill = "#F5FFFA",colour = "#F5FFFA"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25))

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#     Figure 4A NO.5(left to right)
#------------------------------------------------------------

p5<-ggplot(data, aes(x=SubType,
                     y=`Risk Score`,
                     color=SubType),
           palette = c("#FF8C69","#87CEFF")) + 
  geom_boxplot(width=0.5,
               lwd=0.1,
               position=position_dodge(width=0.5))+
  scale_fill_nejm()+
  ggtitle("Test",)+
  geom_jitter(width = 0.2)+
  labs(x="",
       y="Risk Score")+
  geom_signif(comparisons=list(#c("Basal","Her2"),# Add significance tags
    c("Basal","LumA"),
    c("Basal","LumB"),
    #c("Basal","Normal"),
    #c("Her2","LumA"),
    c("Her2","LumB"),
    #c("Her2","Normal"),
    c("LumA","LumB"),
    c("LumA","Normal"),
    c("LumB","Normal")), # Choose which 2 groups you want to add labels to
    step_increase = 0.1,
    test="wilcox.test", # "t-test, comparing two groups (parametric)" = "t.test", "Wilcoxon signed rank test, comparing two groups (nonparametric)" = "wilcox.test"
    map_signif_level=T  # Label style F is a number, T is a * sign
  )+
  theme_par()+
  stat_compare_means(method = "kruskal.test",
                     size = 7,
                     label = "p",
                     # label.x = 1,
                     label.y = 1.8
  )+ #method denotes the test model used, kruskal.test is a nonparametric test and t.test is a parametric test
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        panel.background = element_rect(fill = "#F5FFFA",colour = "#F5FFFA"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle=25, hjust=1, vjust=1))

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#     Figure 4A NO.8(left to right)
#------------------------------------------------------------

p6<-ggplot(data, aes(x=Age,
                     y=`Risk Score`,
                     color=Age),
           palette = c("#FF8C69","#87CEFF")) + 
  geom_boxplot(width=0.5,
               lwd=0.1,
               position=position_dodge(width=0.5))+
  scale_fill_nejm()+
  ggtitle("Test",)+
  geom_jitter(width = 0.2)+
  labs(x="",
       y="Risk Score")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 7,
                     label.y = 0.4)+ 
  theme_par()+
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        panel.background = element_rect(fill = "#F5FFFA",colour = "#F5FFFA"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25))

#############################################    Entire

data<-read.csv("clinic_entire.csv",row.names = 1,check.names = F)
data$Stage<-substr(data$Stage, 6, nchar(data$Stage))

#------------------------------------------------------------
#     Figure 4A NO.3(left to right)
#------------------------------------------------------------

p7<-ggplot(data, aes(x=Stage,
                     y=`Risk Score`,
                     color=Stage),
           palette = c("#FF8C69","#87CEFF")) + 
  geom_boxplot(width=0.5,
               lwd=0.1,
               position=position_dodge(width=0.5))+
  scale_fill_nejm()+
  ggtitle("Entire",)+
  geom_jitter(width = 0.2)+
  labs(x="",
       y="Risk Score")+
  geom_signif(comparisons=list(
    # c("StageI","StageII"),# Add significance tags
    # c("StageI","StageIII"),
    c("StageI","StageIV")
    # c("StageII","StageIII"),
    # c("StageII","StageIV"),
    # c("StageIII","StageIV")
  ), # Choose which 2 groups you want to add labels to
  step_increase = 0.1,
  test="wilcox.test", # "t-test, comparing two groups (parametric)" = "t.test", "Wilcoxon signed rank test, comparing two groups (nonparametric)" = "wilcox.test"
  map_signif_level=T  # Label style F is a number, T is a * sign
  )+
  theme_par()+
  stat_compare_means(method = "kruskal.test",
                     size = 7,
                     label = "p",
                     # label.x = 1,
                     label.y = 1.2
  )+ #method denotes the test model used, kruskal.test is a nonparametric test and t.test is a parametric test
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        panel.background = element_rect(fill = "#F5FFFA",colour = "#F5FFFA"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25))

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#     Figure 4A NO.6(left to right)
#------------------------------------------------------------

p8<-ggplot(data, aes(x=SubType,
                     y=`Risk Score`,
                     #fill=WGD,
                     color=SubType),
           palette = c("#FF8C69","#87CEFF")) + 
  geom_boxplot(width=0.5,
               lwd=0.1,
               position=position_dodge(width=0.5))+
  scale_fill_nejm()+
  ggtitle("Entire",)+
  geom_jitter(width = 0.2)+
  labs(x="",
       y="Risk Score")+
  geom_signif(comparisons=list(#c("Basal","Her2"),# Add significance tags
    #c("Basal","LumA"),
    c("Basal","LumB"),
    c("Basal","Normal"),
    #c("Her2","LumA"),
    c("Her2","LumB"),
    c("Her2","Normal"),
    c("LumA","LumB"),
    c("LumA","Normal"),
    c("LumB","Normal")), # Choose which 2 groups you want to add labels to
    step_increase = 0.1,
    test="wilcox.test", # "t-test, comparing two groups (parametric)" = "t.test", "Wilcoxon signed rank test, comparing two groups (nonparametric)" = "wilcox.test"
    map_signif_level=T  # Label style F is a number, T is a * sign
  )+
  theme_par()+
  stat_compare_means(method = "kruskal.test",
                     size = 7,
                     label = "p",
                     # label.x = 1,
                     label.y = 3
  )+ #method denotes the test model used, kruskal.test is a nonparametric test and t.test is a parametric test
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        panel.background = element_rect(fill = "#F5FFFA",colour = "#F5FFFA"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle=25, hjust=1, vjust=1))

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#     Figure 4A NO.9(left to right)
#------------------------------------------------------------

p9<-ggplot(data, aes(x=Age,
                     y=`Risk Score`,
                     color=Age),
           palette = c("#FF8C69","#87CEFF")) + 
  geom_boxplot(width=0.5,
               lwd=0.1,
               position=position_dodge(width=0.5))+
  scale_fill_nejm()+
  ggtitle("Entire",)+
  geom_jitter(width = 0.2)+
  labs(x="",
       y="Risk Score")+
  stat_compare_means(method = "wilcox.test",
                     label = "p",
                     size = 7,
                     label.y = 0.8)+ 
  theme_par()+
  theme(axis.line =  element_line(colour = "black"),
        plot.title = element_text(hjust = 0,size = 25),
        strip.background = element_rect(color = "white"),
        panel.background = element_rect(fill = "#F5FFFA",colour = "#F5FFFA"),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        axis.text = element_text(size = 25),
        axis.title = element_text(size = 25))

#-----------------------------------------------------------
#-----------------------------------------------------------

#################   Sankey


#############################################    Training
data<-read.csv("clinic_train.csv",row.names = 1, check.names = F)
data$OS[data$OS %in% 0]<-"Alive"
data$OS[data$OS %in% 1]<-"Dead"
data$Group[data$Group %in% "High"]<-"High Group"
data$Group[data$Group %in% "Low"]<-"Low Group "

#Frequency counts in groups
BRCAData <- group_by(data,Age,Group,SubType,OS) %>% summarise(., count = n())

#------------------------------------------------------------
#    Figure 4B NO.1
#------------------------------------------------------------
p10<-ggplot(as.data.frame(BRCAData),
            aes(axis1 = Group, 
                axis2 = OS, 
                axis3 = Age,
                axis4 = SubType,
                y= count))+
  scale_x_discrete(limits = c("Group","OS","Age","SubType"),
                   expand = c(0,.6))+
  geom_alluvium(aes(fill = Group))+
  geom_stratum(width = 1/3, fill = "#DDA0DD", color = "black") +
  theme_minimal() +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)),size = 7) +
  ylab("")+
  scale_fill_manual(values = c("#EE2C2C","#27408B"))+#Custom colors
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=20),
        title = element_text(size=20),
        axis.text.y = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "inches")
  )+
  ggtitle("Training")


#############################################    Test

data<-read.csv("clinic_test.csv",row.names = 1,check.names = F)
data$OS[data$OS %in% 0]<-"Alive"
data$OS[data$OS %in% 1]<-"Dead"
data$Group[data$Group %in% "High"]<-"High Group"
data$Group[data$Group %in% "Low"]<-"Low Group "

#Frequency counts in groups
BRCAData <- group_by(data,Age,Group,SubType,OS) %>% summarise(., count = n())

#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Figure 4A NO.2
#------------------------------------------------------------
p11<-ggplot(as.data.frame(BRCAData),
            aes(axis1 = Group, 
                axis2 = OS, 
                axis3 = Age,
                axis4 = SubType,
                y= count))+
  scale_x_discrete(limits = c("Group","OS","Age","SubType"),
                   expand = c(0,.6))+
  geom_alluvium(aes(fill = Group))+
  geom_stratum(width = 1/3, fill = "#DDA0DD", color = "black") +
  theme_minimal() +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)),size = 7) +
  ylab("")+
  scale_fill_manual(values = c("#EE2C2C","#27408B"))+#Custom colors
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        title = element_text(size=20),
        axis.text.y = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "inches")
  )+
  ggtitle("Test")



#############################################    Entire

data<-read.csv("clinic_entire.csv",row.names = 1,check.names = F)
data$OS[data$OS %in% 0]<-"Alive"
data$OS[data$OS %in% 1]<-"Dead"
data$Group[data$Group %in% "High"]<-"High Group"
data$Group[data$Group %in% "Low"]<-"Low Group "

#Frequency counts in groups
BRCAData <- group_by(data,Age,Group,SubType,OS) %>% summarise(., count = n())
#-----------------------------------------------------------
#-----------------------------------------------------------

#------------------------------------------------------------
#    Figure 4A NO.3
#------------------------------------------------------------
p12<-ggplot(as.data.frame(BRCAData),
            aes(axis1 = Group, 
                axis2 = OS, 
                axis3 = Age,
                axis4 = SubType,
                y= count))+
  scale_x_discrete(limits = c("Group","OS","Age","SubType"),
                   expand = c(0,.6))+
  geom_alluvium(aes(fill = Group))+
  geom_stratum(width = 1/3, fill = "#DDA0DD", color = "black") +
  theme_minimal() +
  geom_label(stat = "stratum", aes(label = after_stat(stratum)),size = 7) +
  ylab("")+
  scale_fill_manual(values = c("#EE2C2C","#27408B"))+#Custom colors
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        title = element_text(size=20),
        axis.text.y = element_blank(),
        legend.text = element_text(size=20),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "inches"),
  )+
  ggtitle("Entire")

#-----------------------------------------------------------
#-----------------------------------------------------------
