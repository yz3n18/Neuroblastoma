#### Figure 4 ####
library("ggpubr")
library('GSVA')
#### Figure 4a/e ####
#### MYCN expression ####

my_comparisons_train <- list(c('MYCN-AMP','Subgroup2'),
                             c('Subgroup1','Subgroup2'),
                             c('Subgroup1','Subgroup3'),
                             c('Subgroup2','Subgroup3'))
merge_data_batch_MYCN_amp_train_train1<-data.frame(merge_data_batch_MYCN_amp,
                                                   train_1_sil,train_2_sil,train_3_sil)
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['MYCN',]
merge_data_batch_MYCN_amp_train_train1<-t(merge_data_batch_MYCN_amp_train_train1)
merge_data_batch_MYCN_amp_train_train1<-as.data.frame(merge_data_batch_MYCN_amp_train_train1)
merge_data_batch_MYCN_amp_train_train1$Group<-'NA'
merge_data_batch_MYCN_amp_train_train1[1:714,2]<-'MYCN-AMP'
merge_data_batch_MYCN_amp_train_train1[715:1748,2]<-'Subgroup1'
merge_data_batch_MYCN_amp_train_train1[1749:2064,2]<-'Subgroup2'
merge_data_batch_MYCN_amp_train_train1[2065:2434,2]<-'Subgroup3'
colnames(merge_data_batch_MYCN_amp_train_train1)<-c('MYCN','Group')

ggviolin(merge_data_batch_MYCN_amp_train_train1, x = "Group", y = 'MYCN',
         fill = "Group",
         add = "boxplot", add.params = list(fill = "white"),palette = c("#D22C6C", "#6688AB", "#97C17E",'#BA86B5'))+
  stat_compare_means(comparisons = my_comparisons_train,size=10)+
  ylab("MYCN expression\n") +
  xlab('')+
  #scale_y_continuous(breaks=seq(0, 1, 0.1))+
  #coord_cartesian(ylim = c(-2500,3000))+
  #geom_hline(yintercept = median(gsva_MYC_matrix[Purity_Test$Group=='MYCN_AMP',]$ImmuneScore), 
  #          size=2,linetype = 2)+
  theme_classic()+
  theme(legend.position = "none",
        
        axis.text.x = element_text(  size=30,color = 'black',hjust = 0.5,vjust=0),
        axis.text.y = element_text( size=30,color = 'black'),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(color="black", size=40,vjust=0),
        plot.margin = unit(c(1,1,1,1),"cm"))
ggsave('Subgroup_AMP_MYCN_expression.pdf',width = 16,height = 12)


#### Figure 4b MYCN score ####
MYCN_signature<-read.csv('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/MYCN-signature.csv',header = T,sep = ',')
MYCN_signature_po<-MYCN_signature[MYCN_signature$Correlation>0,]
# MYCN_signature<-list(MYCN_signature$Gene)
gsva_MYCNscore_matrix<-as.data.frame(matrix(data=NA,
                                      nrow=3084,ncol=2))
colnames(gsva_MYCNscore_matrix)<-c('Sample','MYCN_score_PNAS')
gsva_MYCNscore_matrix[1:714,2]<-t(gsva(as.matrix(merge_data_batch_MYCN_amp),MYCN_signature_po , mx.diff=1,kcdf="Gaussian",
                                 method= "ssgsea",ssgsea.norm=F))

gsva_MYCNscore_matrix[715:1748,2]<-t(gsva(as.matrix(train_1_sil),MYCN_signature_po , mx.diff=1,kcdf="Gaussian",
                                    method= "ssgsea",ssgsea.norm=F))
gsva_MYCNscore_matrix[1749:2126,2]<-t(gsva(as.matrix(test_1_sil),MYCN_signature_po , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYCNscore_matrix[2127:2442,2]<-t(gsva(as.matrix(train_2_sil),MYCN_signature_po , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYCNscore_matrix[2443:2550,2]<-t(gsva(as.matrix(test_2_sil),MYCN_signature_po , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYCNscore_matrix[2551:2920,2]<-t(gsva(as.matrix(train_3_sil),MYCN_signature_po , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYCNscore_matrix[2921:3084,2]<-t(gsva(as.matrix(test_3_sil),MYCN_signature_po , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))

#### Figure 4d ####
#### MYCN correlated genes in Subgroup2 ####
pv<-array(1,c(length(train_2_sil[,1]),1)) #  make a empty list for p-value, 
cor<-array(1,c(length(train_2_sil[,1]),1)) # make a empty list for correlation value
for (i in 1:length(train_2_sil[,1])) #change 
{
  VL<-cor.test(as.numeric(train_2_sil[i,]),
               as.numeric(gsva_MYCNscore_matrix[2127:2442,2]),method = "pearson") #### 
  pv[i,1]<-VL$p.value    #### extract p-value
  cor[i,1]<-VL$estimate #### extract correlation value
}#

fdr<-p.adjust(pv, method = 'BH', n = length(pv)) ## adjust p-value
fdr<-as.data.frame(fdr)
Train_MYCN_correlation_Sub2<-data.frame(rownames(Train_MYCN_correlation_Sub2),
                                      pv,fdr,cor)#,log2_cc_ee) 

####Overlap between DEG of Subgroup2 vs (1+3) and MYCN-pos-correlated genes in Subgroup2
subgroup2DEG<-read.table('GEO+Array/Consensus/Kmean/DEGs.txt',sep = '\t',header = T)
subgroup2DEG<-subgroup2DEG[,c(1,7,11)]
subgroup2DEG_pos<-subset(subgroup2DEG,adj.P.Valsub2< 0.05 & logFCsub2>1)
Train_MYCN_correlation_Sub2<-subset(Train_MYCN_correlation_Sub2,
                                              fdr < 0.05 & cor > 0.3)
Subgroup2_Up_Pos<-intersect(subgroup2DEG_pos$Gene,rownames(Train_MYCN_correlation_Sub2))


