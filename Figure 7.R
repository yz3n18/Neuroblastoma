#### Figure 7 ####
#### GDSC correlation ####
GDSC_AUC<-read.csv('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/GDSC_NB_IC_Wed Sep 16 10_13_57 2020.csv',sep = ',',header = T)
GDSC_AUC$Cosmic.sample.Id<-paste0('DATA.',GDSC_AUC$Cosmic.sample.Id)

GDSC_mRNA_NB<-read.xlsx('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/GDSC_Cell_line_RMA_proc_basalExp_NB.xlsx')
GDSC_mRNA_NB<-GDSC_mRNA_NB[,c(1,3:25)]
colnames(GDSC_mRNA_NB)[1]<-'Gene'
GDSC_mRNA_NB<-na.omit(GDSC_mRNA_NB)


GDSC_mRNA_NB<-aggregate(x=GDSC_mRNA_NB[,2:(length(colnames(GDSC_mRNA_NB)))],by=list(GDSC_mRNA_NB$Gene),FUN=median)
rownames(GDSC_mRNA_NB)<-GDSC_mRNA_NB$Group.1
GDSC_mRNA_NB<-GDSC_mRNA_NB[,-1]

GDSC_AUC_matrix<-as.data.frame(matrix(nrow = 175,ncol =23,data = NA ))
colnames(GDSC_AUC_matrix)<-colnames(GDSC_mRNA_NB)

rownames(GDSC_AUC_matrix)<-unique(GDSC_AUC$Drug.name)
GDSC_AUC_matrix$Drug<-unique(GDSC_AUC$Drug.name)
for (i in 1:length(GDSC_AUC[,1])){
  GDSC_AUC_matrix[GDSC_AUC[i,1],GDSC_AUC[i,4]]<-GDSC_AUC[i,9]
}

GDSC_AUC_matrix[1,1]<-1


GDSC_AUC_mRNA_correlation<-as.data.frame(matrix(nrow = 17417,ncol =525,data = NA )) 
rownames(GDSC_AUC_mRNA_correlation)<-rownames(GDSC_mRNA_NB)
colnames(GDSC_AUC_mRNA_correlation)<-unique(GDSC_AUC$Drug.name)

for (i in 1: length(GDSC_AUC_matrix[,1])){
  gdsc_test<-rbind(GDSC_AUC_matrix[i,1:23],GDSC_mRNA_NB)
  gdsc_test<-as.data.frame(t(na.omit(t(gdsc_test))))
  for (j in 2: length(gdsc_test[,1])){
    VL<-cor.test(as.numeric(gdsc_test[1,]),as.numeric(gdsc_test[j,]),method = 'spearman')
    GDSC_AUC_mRNA_correlation[j-1,(3*i-2)]<-VL$estimate   #### extract p-value
    GDSC_AUC_mRNA_correlation[j-1,(3*i-1)]<-VL$p.value  #### extract correlation value
  }
  GDSC_AUC_mRNA_correlation[,(3*i)]<-p.adjust(GDSC_AUC_mRNA_correlation[,(3*i-1)], 
                                              method = 'BH', 
                                              n = length(GDSC_AUC_mRNA_correlation[,(3*i-2)])) ## adjust p-value
  colnames(GDSC_AUC_mRNA_correlation)[(i*3-2):(i*3)]<-rownames(GDSC_AUC_matrix)[i]
}
write.table(GDSC_AUC_mRNA_correlation,'DSC_AUC_mRNA_correlation_BH_new.txt',sep = '\t')
GDSC_AUC_mRNA_correlation_FDR<-GDSC_AUC_mRNA_correlation

FUN1<-function(x){
  x<-as.numeric(x)
  if (x > 0.01){
    x<-NA
  }
  return(x)
}
for (i in seq(2,525,by=3)){
  for (j in 1:17417){
    GDSC_AUC_mRNA_correlation_FDR[j,i]<-FUN1(GDSC_AUC_mRNA_correlation_FDR[j,i])
  }
  
}
write.table(GDSC_AUC_mRNA_correlation_FDR,'DSC_AUC_mRNA_correlation_FDR_new.txt',sep = '\t')
GDSC_AUC_mRNA_correlation_FDR<-read.table('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/DSC_AUC_mRNA_correlation_FDR_new.txt',sep = '\t',header = T)
GDSC_AUC_mRNA_correlation_FDR_significant<-as.data.frame(matrix(ncol = 5,
                                                                nrow = 20000,data = NA))
GDSC_test_length_start<-1
colnames(GDSC_AUC_mRNA_correlation_FDR)

for (i in seq(1,525,by=3)){
  GDSC_test<-na.omit(GDSC_AUC_mRNA_correlation_FDR[,i:(i+2)])
  GDSC_test_length<-length(GDSC_test[,1])
  GDSC_AUC_mRNA_correlation_FDR_significant[GDSC_test_length_start:(GDSC_test_length+GDSC_test_length_start-1),1]<-rownames(GDSC_test)
  GDSC_AUC_mRNA_correlation_FDR_significant[GDSC_test_length_start:(GDSC_test_length+GDSC_test_length_start-1),2:4]<-GDSC_test
  GDSC_AUC_mRNA_correlation_FDR_significant[GDSC_test_length_start:(GDSC_test_length+GDSC_test_length_start-1),5]<-rep(colnames(GDSC_AUC_mRNA_correlation_FDR)[i],
                                                                                                                       GDSC_test_length)
  GDSC_test_length_start<-GDSC_test_length+GDSC_test_length_start
  
  
}
write.csv(GDSC_AUC_mRNA_correlation_FDR_significant,'GDSC_AUC_mRNA_correlation_FDR_significant_new.csv')
rm(GDSC_mRNA,GDSC_AUC,GDSC_AUC_matrix,GDSC_AUC_mRNA_correlation,GDSC_AUC_mRNA_correlation_FDR,GDSC_mRNA_NB,
   GDSC_test,GDSC_test_length,GDSC_test_length_start,gdsc_test)

#### Overlap with DEG ####
DEGs<-read.table('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/DEGs.txt',sep = '\t',header = T)
DEGs<-DEGs[,c(1,2,6,7,11,12,16)]

DEGs_Sub1<-subset(DEGs,adj.P.Valsub1< 0.05 & (logFCsub1>0.5849625 ))
DEGs_Sub2<-subset(DEGs,adj.P.Valsub2< 0.05 & (logFCsub2>0.5849625 ))
DEGs_Sub3<-subset(DEGs,adj.P.Valsub3< 0.05 & (logFCsub3>0.5849625))
GDSCtest<-GDSC_AUC_mRNA_correlation_FDR_significant
colnames(GDSCtest)<-c('Gene','Spearman Correlation Coefficent','P-value','FDR','Drug')
GDSCtest<-merge(GDSCtest,DEGs_Sub1,by='Gene')
write.table(GDSCtest,'GDSCtest_subgroup1.txt',sep = '\t') ## GDSCdrug_DEG_correlation_subgroup.txt

GDSCtest<-GDSC_AUC_mRNA_correlation_FDR_significant
colnames(GDSCtest)<-c('Gene','Spearman Correlation Coefficent','P-value','FDR','Drug')
GDSCtest<-merge(GDSCtest,DEGs_Sub2,by='Gene')
write.table(GDSCtest,'GDSCtest_subgroup2.txt',sep = '\t') ##GDSCdrug_DEG_correlation_subgroup.txt

GDSCtest<-GDSC_AUC_mRNA_correlation_FDR_significant
colnames(GDSCtest)<-c('Gene','Spearman Correlation Coefficent','P-value','FDR','Drug')
GDSCtest<-merge(GDSCtest,IPA_DEGs_Sub3,by='Gene')
write.table(GDSCtest,'GDSCtest_subgroup3.txt',sep = '\t') ## GDSCdrug_DEG_correlation_subgroup.txt

#### Figure 7a/b ####

devtools::install_github("zzwch/crosslinks") 
install.packages("aplot")
install.packages("dplyr")
library(ggplot2)
library(aplot)
library(dplyr)
library(crosslinks)
sub1Data <- read.table("/Volumes/My Passport/2021/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/Figure folder/Figure7 Drug/Figure7-2021/Drug_Subgroup3_sensitive/easy_input_subgroup1.txt", header = T,sep = '\t')
sub2Data <- read.table("/Volumes/My Passport/2021/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/Figure folder/Figure7 Drug/Figure7-2021/Drug_Subgroup3_sensitive/easy_input_subgroup2.txt", header = T)
sub3Data <- read.table("/Volumes/My Passport/2021/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/Figure folder/Figure7 Drug/Figure7-2021/Drug_Subgroup3_sensitive/easy_input_subgroup3.txt", header = T)
nodes <- read.table("/Volumes/My Passport/2021/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/Figure folder/Figure7 Drug/Figure7-2021/Drug_Subgroup3_sensitive/easy_input_nodes.txt", header = T,
                    sep='\t')
# Edge between each dots
edges <- read.table("/Volumes/My Passport/2021/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/Figure folder/Figure7 Drug/Figure7-2021/Drug_Subgroup3_sensitive/easy_input_edges.txt", header = T,
                    sep = '\t',stringsAsFactors = F)
top <- sub2Data$sub2gene[order(sub2Data$lfc, decreasing = T)] # set the order according to log fold change
bottom <- sub3Data$sub3gene[order(sub3Data$lfc, decreasing = T)] # set the order according to log fold change
left <- sub1Data$Sub1gene[order(sub1Data$lfc, decreasing = T)] #  set the order according to log fold change
right <- nodes$id[nodes$type == "Drug"] # 
columns <- list(
  Drug = nodes$id[nodes$type == "Drug"], 
  Target = nodes$id[nodes$type == "Tar"], 
  Pathway = nodes$id[nodes$type == "Path"]
)

columns$Drug<-factor(columns$Drug,levels = columns$Drug)
columns$Target<-factor(columns$Target,levels = columns$Target)
columns$Pathway<-factor(columns$Pathway,levels = columns$Pathway)

height <- 8
width <- 12
flank_mult <- c(top = 0.1, bottom = 0.1,left = 0.1, right = 0.1)

# Three bar plot
topAnn <- sub2Data %>%
  mutate(gene_f = factor(sub2gene, top)) %>%
  ggplot(mapping = aes(
    newCoord(gene_f, width, flank_mult[['top']]),
    lfc)) +
  geom_bar(fill = "#00BFC4",
           stat = "identity",
           width = 0.5) + # 
  geom_text(mapping = aes(label = sub2gene), angle = 90, hjust =-0.2, fontface ="italic") +
  labs(x = NULL, y = "log2(Fold Change)") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"))
topAnn

# Bottom plot
botAnn <- sub3Data %>%
  mutate(meth_f = factor(sub3gene, bottom)) %>%
  ggplot(mapping = aes(newCoord(meth_f, width, flank_mult[['bottom']]),
                       lfc)) +
  geom_bar(fill = "#C77CFF",
           stat = "identity",
           width = 0.5) +
  geom_text(mapping = aes(label = sub3gene), angle = 90, hjust =0, fontface ="italic") +
  labs(x = NULL, y = "Difference") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"))+scale_y_reverse()
botAnn

# left plot
lftAnn <- sub1Data %>%
  mutate(mir_f = factor(Sub1gene, left)) %>%
  ggplot(mapping = aes(newCoord(mir_f, height, flank_mult[['left']]),
                       lfc)) +
  geom_bar(fill = "#7CAE00",
           stat = "identity", 
           width = 0.5) +
  geom_text(mapping = aes(y = 0,label = Sub1gene), 
            #fontface ="italic", 
            vjust = -1, 
            hjust = 1) +
  labs(x = NULL, y = "log2(Fold Change)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,1,1,1), "cm"))+
  coord_flip()+
  scale_y_reverse()

lftAnn
squareCross(edges, top ,bottom, left, right, height = 8, width = 12, flank_mult = 0.1)%>%
  insert_top(topAnn, 0.4) %>%
  insert_bottom(botAnn, 0.4) %>%
  insert_left(lftAnn, 0.4)
# Merge line with bar plots
p <- (squareCross(edges, top ,bottom, left, right, height = 8, width = 12, flank_mult = 0.1) +
        scale_color_manual(values = c("darkseagreen3",'goldenrod1','goldenrod1'))) %>%
  insert_top(topAnn, 0.4) %>%
  insert_bottom(botAnn, 0.4) %>%
  insert_left(lftAnn, 0.4) %>%
  
  # 
  insert_right(columnCross(edges, columns, height = 8, flank_mult = 0.1, segment_shrink = 0.2) +
                 scale_color_manual(values = c("red", "red")) , 0.6)
ggsave(plot = p, filename = "Figure7a_sub2.pdf", width = 20, height = 15)
