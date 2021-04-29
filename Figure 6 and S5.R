#### Figure 6 ####
#### Immune pathway in Subgroup3 ####
Immune_pathway<-read.csv('The immune related pathway and genes.csv',header = T,sep = ',')
dim(Immune_pathway)

gsva_Immunepathway_matrix<-as.data.frame(matrix(data=NA,
                                                nrow=3084,ncol=49))
colnames(gsva_Immunepathway_matrix)[1:46]<-colnames(Immune_pathway)
colnames(gsva_Immunepathway_matrix)[47:49]<-c('Group','Group-MYCN','Group-subgroup')

i<-'CDK4'
c(mean(train_1[i,]),mean(train_2[i,]),mean(train_3[i,]))

for (i in 1:46){
  gsva_Immunepathway_matrix[1:714,i]<-t(gsva(as.matrix(merge_data_batch_MYCN_amp),list(Immune_pathway[,i]) , mx.diff=1,kcdf="Gaussian",
                                             method= "ssgsea",ssgsea.norm=F))
  gsva_Immunepathway_matrix[715:1748,i]<-t(gsva(as.matrix(train_1),list(Immune_pathway[,i]) , mx.diff=1,kcdf="Gaussian",
                                                method= "ssgsea",ssgsea.norm=F))
  gsva_Immunepathway_matrix[1749:2126,i]<-t(gsva(as.matrix(test_1),list(Immune_pathway[,i]) , mx.diff=1,kcdf="Gaussian",
                                                 method= "ssgsea",ssgsea.norm=F))
  gsva_Immunepathway_matrix[2127:2442,i]<-t(gsva(as.matrix(train_2),list(Immune_pathway[,i]) , mx.diff=1,kcdf="Gaussian",
                                                 method= "ssgsea",ssgsea.norm=F))
  gsva_Immunepathway_matrix[2443:2550,i]<-t(gsva(as.matrix(test_2),list(Immune_pathway[,i]) , mx.diff=1,kcdf="Gaussian",
                                                 method= "ssgsea",ssgsea.norm=F))
  gsva_Immunepathway_matrix[2551:2920,i]<-t(gsva(as.matrix(train_3),list(Immune_pathway[,i]) , mx.diff=1,kcdf="Gaussian",
                                                 method= "ssgsea",ssgsea.norm=F))
  gsva_Immunepathway_matrix[2921:3084,i]<-t(gsva(as.matrix(test_3),list(Immune_pathway[,i]) , mx.diff=1,kcdf="Gaussian",
                                                 method= "ssgsea",ssgsea.norm=F))
  
}
gsva_Immunepathway_matrix[1:714,47]<-'Neuroblastoma'
gsva_Immunepathway_matrix[1:714,48]<-'MYCN-AMP'
gsva_Immunepathway_matrix[1:714,49]<-'MYCN-AMP'
gsva_Immunepathway_matrix[715:2126,47]<-'Neuroblastoma'
gsva_Immunepathway_matrix[715:2126,48]<-'MYCN-Normal'
gsva_Immunepathway_matrix[715:2126,49]<-'Subgroup1'
gsva_Immunepathway_matrix[2127:2550,47]<-'Neuroblastoma'
gsva_Immunepathway_matrix[2127:2550,48]<-'MYCN-Normal'
gsva_Immunepathway_matrix[2127:2550,49]<-'Subgroup2'
gsva_Immunepathway_matrix[2551:3084,47]<-'Neuroblastoma'
gsva_Immunepathway_matrix[2551:3084,48]<-'MYCN-Normal'
gsva_Immunepathway_matrix[2551:3084,49]<-'Subgroup3'

gsva_Immunepathway_matrix[1:714,50]<-'MYCN-AMP'
gsva_Immunepathway_matrix[715:1748,50]<-'Train'
gsva_Immunepathway_matrix[1749:2126,50]<-'Test'
gsva_Immunepathway_matrix[2127:2442,50]<-'Train'
gsva_Immunepathway_matrix[2443:2550,50]<-'Test'
gsva_Immunepathway_matrix[2551:2920,50]<-'Train'
gsva_Immunepathway_matrix[2921:3084,50]<-'Test'

gsva_Immunepathway_matrix_train<-gsva_Immunepathway_matrix[gsva_Immunepathway_matrix$Train_Test!='Test',]
gsva_Immunepathway_matrix_train<-gsva_Immunepathway_matrix_train[,1:46]
gsva_Immunepathway_matrix_train<-t(gsva_Immunepathway_matrix_train)
n=t(scale(t(gsva_Immunepathway_matrix_train)))
n[n>3]=3
n[n< -3]= -3

ac<-data.frame(c(rep('MYCN-AMP',714),rep('Subgroup1',1034),rep('Subgroup2',316), # AMP+Train
                 rep('Subgroup3',370)))

colnames(ac)<-'Molecular subtype'
#ac$`Molecular subtype` <- factor(ac$`Molecular subtype`, levels = c("MYCN-AMP", "Subgroup1",'Subgroup2','Subgroup3'))
rownames(ac)=colnames(n)
ac$Blank<-'1'
ac<-as.matrix(t(ac))
rownames(ac)=colnames(n_order)


pdf("heatmap_Immunepathway_order.pdf", width = 15, height = 10)
ha = HeatmapAnnotation(`Molecular subtype` =ac[1,],
                       col = list(`Molecular subtype` = c('MYCN-AMP' = "#F8766D", 'Subgroup1' = "#7CAE00",
                                                          'Subgroup2'='#00BFC4','Subgroup3'='#C77CFF')),
                       annotation_legend_param = list(
                         `Molecular subtype` = list(title = "Molecular subtype")))
draw(Heatmap(as.matrix(n), name = "Immune pathway", 
             column_km = 1,column_km_repeats = 1000,
             clustering_distance_columns = "euclidean",
             cluster_columns = FALSE,
             show_row_names = T,show_column_names  = FALSE,
             top_annotation = ha,col=colorbar,show_heatmap_legend = T,
             show_parent_dend_line = FALSE),padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()


#### Figure 6b/c ####
merge_data_batch_MYCN_amp_train_train1<-data.frame(merge_data_batch_MYCN_amp,train_1_sil,train_2_sil,train_3_sil)

MHC_CYT<-data.frame(matrix(data=NA,ncol=4,nrow=2434))
colnames(MHC_CYT)<-c('Sample','Group','MHC-1','CYT')
MHC_CYT$Sample<-colnames(merge_data_batch_MYCN_amp_train_train1)
MHC_CYT[1:714,2]<-'MYCN-AMP'
MHC_CYT[715:1748,2]<-'Subgroup1'
MHC_CYT[1749:2064,2]<-'Subgroup2'
MHC_CYT[2065:2434,2]<-'Subgroup3'

merge_data_batch_MYCN_amp_train_train1_MHC<-merge_data_batch_MYCN_amp_train_train1[c('HLA-A', 'HLA-B', 
                                                                                     'HLA-C', 'TAP1', 'TAP2', 
                                                                                     'NLRC5', 'PSMB9',
                                                                                     'PSMB8','B2M'),]

MHC_CYT$`MHC-1`<-apply(merge_data_batch_MYCN_amp_train_train1_MHC[1:9,],2,mean)
#### Cytolytic activity was calculated as the geometric mean of GZMA and PRF1 
merge_data_batch_MYCN_amp_train_train1_MHC<-merge_data_batch_MYCN_amp_train_train1[c('GZMA','PRF1'),]
FUN1<-function(x){
  sqrt(x[1]*x[2])
}
MHC_CYT$CYT<-apply(merge_data_batch_MYCN_amp_train_train1_MHC,2,FUN1)
merge_data_batch_MYCN_amp_train_train1_MHC[,1]

ggplot(MHC_CYT, aes(CYT,colour =Group ),size=5) + stat_ecdf(geom = "step")+
  labs(
    y = "Cumulative Density Function (CDF)\n",
    x='\nCYT score')+
  theme_classic()+
  theme(axis.text.x = element_text(  size=25,color = 'black'),
        axis.text.y = element_text(  size=25,color = 'black'),
        axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.title.y = element_text(color="black", size=25),
        axis.title.x = element_text(color="black", size=25),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.text = element_text( size = 20),
        legend.title=element_blank())
ggsave('CYT_AMP_Train.pdf',width = 10,height = 8)

ggplot(MHC_CYT, aes(`MHC-1`,colour =Group ),size=5) + stat_ecdf(geom = "step")+
  labs(
    y = "Cumulative Density Function (CDF)\n",
    x='\nMHC-1 score')+
  theme_classic()+
  theme(axis.text.x = element_text(  size=25,color = 'black'),
        axis.text.y = element_text(  size=25,color = 'black'),
        axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
        axis.title.y = element_text(color="black", size=25),
        axis.title.x = element_text(color="black", size=25),
        plot.margin = unit(c(1,1,1,1), "cm"),
        legend.text = element_text( size = 20),
        legend.title=element_blank())
ggsave('MHC-1_AMP_Train.pdf',width = 10,height = 8)

#### Kolmogorov-Smirnov Tests ####
library(dgof)
ks_1<-dgof::ks.test(MHC_CYT[1:714,3], MHC_CYT[2065:2434,3])
ks.test(MHC_CYT[715:1748,3], MHC_CYT[2065:2434,3])
ks.test(MHC_CYT[1749:2064,3], MHC_CYT[2065:2434,3])

ks.test(MHC_CYT[1:714,4], MHC_CYT[2065:2434,4])
ks.test(MHC_CYT[715:1748,4], MHC_CYT[2065:2434,4])
ks.test(MHC_CYT[1749:2064,4], MHC_CYT[2065:2434,4])
#### Figure 6d & Supplementary Figure 5a/b ####
#### Tumor purity ####
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(stringr)

table(colnames(merge_data_batch))
FUN1<-function(x){
  x<-strsplit(x,'_')[[1]][1]
}
a<-sapply(colnames(merge_data_batch),FUN1)
dim(merge_data_batch)
table(a)
write.table(merge_data_batch[,1:393],'GSE120572_RMBATCH.txt',sep = '\t') #GSE120572 agilent 393
write.table(merge_data_batch[,393:445],'GSE112447_RMBATCH.txt',sep = '\t')#GSE112447 Agilent 52
write.table(merge_data_batch[,446:550],'GSE73517_RMBATCH.txt',sep = '\t')#GSE73517 Agilent 105
write.table(merge_data_batch[,551:1043],'GSE49710_RMBATCH.txt',sep = '\t')#GSE49710 Agilent 493
write.table(merge_data_batch[,1044:1686],'GSE45547_RMBATCH.txt',sep = '\t')#GSE45547 Agilent 643
write.table(merge_data_batch[,1687:1786],'GSE19274_RMBATCH.txt',sep = '\t')#GSE19274 Illumina 100
write.table(merge_data_batch[,1787:1837],'GSE16237_RMBATCH.txt',sep = '\t')#GSE16237 Affymetrix 51
write.table(merge_data_batch[,1838:1866],'GSE14880_RMBATCH.txt',sep = '\t')#GSE14880 Affymetrix 29
write.table(merge_data_batch[,1867:1896],'GSE13136_RMBATCH.txt',sep = '\t')#GSE13136 Affymetrix 30
write.table(merge_data_batch[,1897:1943],'GSE27608_RMBATCH.txt',sep = '\t')#GSE27608 Affymetrix 47
write.table(merge_data_batch[,1944:1983],'GSE21713_RMBATCH.txt',sep = '\t')#GSE21713 Affymetrix 40
write.table(merge_data_batch[,1984:2260],'GSE85047_RMBATCH.txt',sep = '\t')#GSE85047 Affymetrix 277
write.table(merge_data_batch[,2261:2348],'GSE16476_RMBATCH.txt',sep = '\t')#GSE16476 Affymetrix 88
write.table(merge_data_batch[,2349:2375],'GSE79910_RMBATCH.txt',sep = '\t')#GSE79910 Agilent 27
write.table(merge_data_batch[,2376:3058],'EMTAB1781_RMBATCH.txt',sep = '\t')#E_MTAB_1781 Agilent 683
write.table(merge_data_batch[,3059:3280],'EMTAB8248_RMBATCH.txt',sep = '\t')#E_MTAB_8248 Agilent 222
write.table(merge_data_batch[,3281:3757],'EMTAB179_RMBATCH.txt',sep = '\t')#E_MTAB_179 Agilent 477
write.table(merge_data_batch[,3758:3799],'EMEXP3517_RMBATCH.txt',sep = '\t')#EMEXP3517 Affymetrix 42
dim(merge_data_batch)
filterCommonGenes(
  input.f = 'EMEXP3517_RMBATCH.txt',
  output.f = "OV_10412genes_10.gct",
  id="GeneSymbol"
)
estimateScore(
  "OV_10412genes_10.gct",
  "EMEXP3517_RMBATCH_Affymetrix.gct",
  platform="affymetrix"
)

Purity_All<-read.csv('/Users/lefan/Rcode/Neuroblastoma/Hu/GEO+Array/Tumorpurity.csv',header = T,sep = ',')
rownames(Purity_All)<-Purity_All$Sample
Purity_All<-Purity_All[,-5] # for score
Purity_All<-na.omit(Purity_All) # for purity
Purity_All<-Purity_All[c(intersect(rownames(Purity_All),colnames(merge_data_batch_MYCN_normal_train_estimate)),
                         intersect(rownames(Purity_All),colnames(merge_data_batch_MYCN_normal_test_estimate)),
                         intersect(rownames(Purity_All),colnames(merge_data_batch_MYCN_amp))),]
Purity_All$Group<-'MYCN_AMP'
Purity_All[intersect(rownames(Purity_All),colnames(train_1)),5]<-'Subgroup_1' # 6 for purity / 5 for score
Purity_All[intersect(rownames(Purity_All),colnames(train_2)),5]<-'Subgroup_2'
Purity_All[intersect(rownames(Purity_All),colnames(train_3)),5]<-'Subgroup_3'
Purity_All[intersect(rownames(Purity_All),colnames(test_1)),5]<-'Subgroup_1'
Purity_All[intersect(rownames(Purity_All),colnames(test_2)),5]<-'Subgroup_2'
Purity_All[intersect(rownames(Purity_All),colnames(test_3)),5]<-'Subgroup_3'

merge_data_batch_MYCN_normal_train_estimate<-data.frame(train_1_sil,train_2_sil,train_3_sil)

Purity_Train<-Purity_All[c(intersect(rownames(Purity_All),colnames(merge_data_batch_MYCN_normal_train_estimate)),
                           intersect(rownames(Purity_All),colnames(merge_data_batch_MYCN_amp))),]
Purity_Train<-na.omit(Purity_Train)
Purity_Train$Group<-'MYCN_AMP'
Purity_Train[colnames(test_1_sil),6]<-'Subgroup_1'
Purity_Train<-na.omit(Purity_Train)
Purity_Train[colnames(test_2_sil),6]<-'Subgroup_2'
Purity_Train<-na.omit(Purity_Train)
Purity_Train[colnames(test_3_sil),6]<-'Subgroup_3'
Purity_Train<-na.omit(Purity_Train)
dim(Purity_Train)

merge_data_batch_MYCN_normal_test_estimate<-data.frame(test_1_sil,test_2_sil,test_3_sil)
Purity_Test<-Purity_All[c(intersect(rownames(Purity_All),colnames(merge_data_batch_MYCN_normal_test_estimate)),
                          intersect(rownames(Purity_All),colnames(merge_data_batch_MYCN_amp))),]
Purity_Test<-na.omit(Purity_Test)

Purity_Test$Group<-'MYCN_AMP'
Purity_Test[colnames(test_1_sil),6]<-'Subgroup_1'
Purity_Test<-na.omit(Purity_Test)
Purity_Test[colnames(test_2_sil),6]<-'Subgroup_2'
Purity_Test<-na.omit(Purity_Test)
Purity_Test[colnames(test_3_sil),6]<-'Subgroup_3'
Purity_Test<-na.omit(Purity_Test)
rownames(Purity_Train)<-Purity_Train$Sample
library(ggplot2)
library(tidyverse)
library(ggpubr)
my_comparisons_train <- list( c("Subgroup_1", "MYCN_AMP"), 
                              c("Subgroup_2", "MYCN_AMP"),
                              c("Subgroup_3", "MYCN_AMP"))
colnames(Purity_Test)
Purity_Test$Group <- factor(Purity_Test$Group,levels = c("MYCN_AMP", "Subgroup_1", "Subgroup_2",
                                                         "Subgroup_3"))
Purity_Train$Group <- factor(Purity_Train$Group,levels = c("MYCN_AMP", "Subgroup_1", "Subgroup_2",
                                                           "Subgroup_3"))
Purity_All$Group <- factor(Purity_All$Group,levels = c("MYCN_AMP", "Subgroup_1", "Subgroup_2",
                                                       "Subgroup_3"))
a<-ggviolin(Purity_Test, x = "Group", y = "ImmuneScore", fill = "Group",
            add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons_train,size=10)+
  ylab("ImmuneScore\n") +
  xlab('')+
  #scale_y_continuous(breaks=seq(0, 1, 0.1))+
  #coord_cartesian(ylim = c(-2500,3000))+
  geom_hline(yintercept = median(Purity_Test[Purity_Test$Group=='MYCN_AMP',]$ImmuneScore), 
             size=2,linetype = 2)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text( face="bold", size=30,color = 'black',angle=45,hjust = 1),
        axis.text.y = element_text( face="bold", size=30,color = 'black'),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(color="black", size=40, face="bold"))
#ggsave('Train_TumorPurity.pdf',width = 10,height = 8) 

b<-ggviolin(Purity_Train, x = "Group", y = "ImmuneScore", fill = "Group",
            add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons_train,size=10)+
  ylab("Immune Score\n") +
  xlab('')+
  #scale_y_continuous(breaks=seq(0, 1, 0.1))+
  #coord_cartesian(ylim = c(-2500,3000))+
  geom_hline(yintercept = median(Purity_Train[Purity_Train$Group=='MYCN_AMP',]$ImmuneScore), 
             size=2,linetype = 2)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text( face="bold", size=30,color = 'black',angle=45,hjust = 1),
        axis.text.y = element_text( face="bold", size=30,color = 'black'),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(color="black", size=40, face="bold"))
c<-ggviolin(Purity_All, x = "Group", y = "ImmuneScore", fill = "Group",
            add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons_train,size=10)+
  ylab("Immune Score\n") +
  xlab('')+
  #scale_y_continuous(breaks=seq(0, 1, 0.1))+
  #coord_cartesian(ylim = c(-2500,3000))+
  geom_hline(yintercept = median(Purity_All[Purity_All$Group=='MYCN_AMP',]$ImmuneScore), 
             size=2,linetype = 2)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text( face="bold", size=30,color = 'black',angle=45,hjust = 1),
        axis.text.y = element_text( face="bold", size=30,color = 'black'),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(color="black", size=40, face="bold"))
d<-ggarrange(a, b, c,
             ncol = 3, nrow = 1)+
  theme(plot.margin = unit(c(4,2,2,2), "cm"))
ggsave(d,filename = 'Tumor ImmuneScore Train Test_silhouette.pdf',width = 36,height = 12) 


#### Single cell ####
library(dplyr)
library(Seurat)
library(patchwork)
setwd('/home/yz3n18/single-cell/')
options(future.globals.maxSize = 200 * 1024^3)

T174 <- Read10X(data.dir = "./T175/")
T188 <- Read10X(data.dir = "./T188/")
T214 <- Read10X(data.dir = "./T214/")
T174<-CreateSeuratObject(counts = T174, 
                         project = "T174",min.cells = 10, min.features = 200)
T188<-CreateSeuratObject(counts = T188, 
                         project = "T188",min.cells = 10, min.features = 200)
T214<-CreateSeuratObject(counts = T214, 
                         project = "T214",min.cells = 10, min.features = 200)

T174 [["percent.mt"]] <- PercentageFeatureSet(T174 , pattern = "^mt-")
T174 <- subset(T174,subset = nFeature_RNA > 200 &  nCount_RNA > 500 &  nCount_RNA < 20000 & percent.mt < 10)
T188 [["percent.mt"]] <- PercentageFeatureSet(T188 , pattern = "^mt-")
T188 <- subset(T188,subset = nFeature_RNA > 200 &  nCount_RNA > 500 &  nCount_RNA < 20000 & percent.mt < 10)
T214 [["percent.mt"]] <- PercentageFeatureSet(T214 , pattern = "^mt-")
T214 <- subset(T214,subset = nFeature_RNA > 200 &  nCount_RNA > 500 &  nCount_RNA < 20000 & percent.mt < 10)
T174$group<-'V174'
T188$group<-'V188'
T214$group<-'V214'
T174 <- SCTransform(T174, vars.to.regress = "percent.mt", verbose = FALSE,variable.features.n=5000)
T188 <- SCTransform(T188, vars.to.regress = "percent.mt", verbose = FALSE,variable.features.n=5000)
T214 <- SCTransform(T214, vars.to.regress = "percent.mt", verbose = FALSE,variable.features.n=5000)

T174 <- RenameCells(object = T174, add.cell.id = "T174")
T188 <- RenameCells(object = T188, add.cell.id = "T188")
T214 <- RenameCells(object = T214, add.cell.id = "T214")

pancreas_12.list<-list(T174,T188,T214)
pancreas_12.features <- SelectIntegrationFeatures(object.list = pancreas_12.list, nfeatures = 5000)

pancreas_12.list <- PrepSCTIntegration(object.list = pancreas_12.list, 
                                       anchor.features = pancreas_12.features, 
                                       verbose = FALSE)

pancreas_12.anchors <- FindIntegrationAnchors(object.list = pancreas_12.list, normalization.method = "SCT", 
                                              anchor.features = pancreas_12.features, verbose = FALSE)
pancreas_12.integrated <- IntegrateData(anchorset = pancreas_12.anchors, normalization.method = "SCT", 
                                        verbose = FALSE)
set.seed(123)
pancreas_12.integrated <- RunPCA(pancreas_12.integrated, npcs = 40,verbose = FALSE)
set.seed(123)
pancreas_12.integrated <- RunUMAP(pancreas_12.integrated,reduction = "pca", dims = 1:40)
#save(pancreas_12.integrated,file='MYCN-non-pancreas_12integrated_afterUMAP.Rdata')
setwd('/home/yz3n18/single-cell/clusterMYCNnon')
pancreas_12.integrated <- FindNeighbors(pancreas_12.integrated,reduction = "pca", dims = 1:40, verbose = FALSE)
save(pancreas_12.integrated,file='MYCN-non-pancreas_12.integrated_Find.Rdata')
rm(pancreas_12.anchors,pancreas_12.list,pancreas_12.features,T174,T188,T214)
dim.use<-1:40
i<-0.9
pancreas.integrated <- FindClusters(pancreas_12.integrated, verbose = FALSE,resolution =i )
pdf(paste0('1111-Integrated_Resolution_',i,"_CellCluster-UMAPPlot_",max(40),"PC.pdf"),width = 10,height = 8)
DimPlot(pancreas.integrated,reduction = "umap", pt.size=0.5,label = T,label.size = 8)
dev.off()
#### Supplementary Figure 5c ####
pdf(paste0('1111-Integrated_Resolution_',i,"_CellCluster-UMAPPlot_",max(40),"PC-group.pdf"),width = 48,height = 8)
DimPlot(pancreas.integrated, reduction = "umap", split.by = "group")
dev.off()
all.markers <- FindAllMarkers(pancreas.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(all.markers,file=paste0('0311-Integrated_Resolution_',i,"_marker_genes_UMAP_",max(dim.use),"PC.txt"),sep="\t",quote = F)
pdf(paste0('1111-Integrated_Resolution_',i,"_CellCluster-UMAPPlot_",max(40),"PC-Vlnplot.pdf"),width = 80,height = 160)
VlnPlot(pancreas.integrated, features = c('CD2','CD3D','CD3E','CD3G','CD7','CD52',"CD79A", "CD79B",'CD19','MS4A1','CD52','CD14','CD68','AIF1','CSF1R','TYROBP',         'FCER1G','DCN','TAGLN','COL3A1','COL1A1','COL6A1','ENG','CLDN5','VMF','CDH5','IL3RA','CLEC4C','SOX10','S100B','PHOX2B','PHOX2A','TH','SYP','CHGA','CHGB','DBH','GATA2'))
dev.off()

pdf(paste0('1111-Integrated_Resolution_',i,"_CellCluster-UMAPPlot_",max(40),"PC-Featureplot.pdf"),width = 80,height = 160)
FeaturePlot(pancreas.integrated, features = c('CD2','CD3D','CD3E','CD3G','CD7','CD52',"CD79A", "CD79B",'CD19','MS4A1','CD52','CD14','CD68','AIF1','CSF1R','TYROBP',                'FCER1G','DCN','TAGLN','COL3A1','COL1A1','COL6A1','ENG','CLDN5','VMF','CDH5','IL3RA','CLEC4C','SOX10','S100B','PHOX2B','PHOX2A',                'TH','SYP','CHGA','CHGB','DBH','GATA2'))
dev.off()




