#### Merged_data.Rdata contains merge_data and MYCN_conditions
#### merge_data is the merged file of 18 MYCN-normal datasets | MYCN_conditions contains MYCN_status of all samples

load('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Merged_data.Rdata')

#### PCA ####
library(factoextra)
data_t <- t(as.matrix(merge_data))
variableL <- ncol(data_t)
pca <- prcomp(data_t[,1:variableL], scale=T)
fviz_pca_ind(pca, repel=T)  
data_t<-as.data.frame(data_t)
data_t$conditions<-as.character(mycn_state_newe)
pdf('PCA_withoutremovebatch.pdf',width = 8,height = 8)
fviz_pca_ind(pca, col.ind=MYCN_conditions$Dataset, mean.point=F, addEllipses = T, legend.title="Dataset",
             geom.ind = "point")+scale_shape_manual(values=seq(0,17))+theme_classic()
dev.off()


#### remove batch effect 
library(sva)
table(MYCN_conditions$Dataset)

### create batch group

MYCN_conditions$status<-c(rep(1,132),rep(2,52),rep(3,101),
                          rep(4,33),rep(5,94),rep(6,100),
                          rep(7,51),rep(8,30),rep(9,47),
                          rep(10,277),rep(11,88),rep(12,27),
                          rep(13,41),rep(14,184),rep(15,262),rep(16,47))


batch<-MYCN_conditions$status

### remove batch
merge_data_batch<- ComBat(dat=as.matrix(merge_data), batch=batch)

#### Supplementary Figure 3a ####
pdf('Boxplot_without_removebatch.pdf',width = 40,height = 8)
par(mar=c(2,8,2,2))
boxplot(merge_data,names=NA,cex.axis=4)
dev.off()
pdf('Boxplot_removebatch.pdf',width = 40,height = 8)
par(mar=c(2,8,2,2))
boxplot(merge_data_batch,names=NA,cex.axis=4)
dev.off()


FUN1<-function(x){
  if (x == 0){
    x<-'MYCN-normal'
  } else{
    x<-'MYCN-amplification'
  }
}
dim(MYCN_conditions)
MYCN_conditions_new<-MYCN_conditions
MYCN_conditions_new$MYCN_conditions<-sapply(MYCN_conditions$MYCN_conditions, FUN1)
MYCN_conditions$MYCN_conditions

#### Supplementary Figure 1b ####
pdf('PCA_removebatch_MYCNstatus.pdf',width = 8,height = 8)
fviz_pca_ind(pca_batch, col.ind=MYCN_conditions_new$MYCN_conditions, mean.point=F, addEllipses = T, 
             legend.title="Status",
             geom.ind = "point")+scale_shape_manual(values=seq(0,17))+theme_classic()
dev.off()


#### extract MYCN-amp from merge_batch####
merge_data_batch_MYCN_amp<-merge_data_batch[,MYCN_conditions_new$MYCN_conditions=='MYCN-amplification']
MYCN_conditions_new_amp<-MYCN_conditions_new[MYCN_conditions_new$MYCN_conditions=='MYCN-amplification',]

#### extract MYCN-normal from merge_batch ####
merge_data_batch_MYCN_normal<-merge_data_batch[,MYCN_conditions_new$MYCN_conditions=='MYCN-normal']
MYCN_conditions_new_normal<-MYCN_conditions_new[MYCN_conditions_new$MYCN_conditions=='MYCN-normal',]

#### Randomly divide data into train and test ####
{install.packages("caret")
  library(caret)
  set.seed(123)
  ind <- MYCN_conditions_new_normal$MYCN_conditions %>% createDataPartition(p = 0.7, list = FALSE)
  merge_data_batch_MYCN_normal_train<-merge_data_batch_MYCN_normal[,ind]
  merge_data_batch_MYCN_normal_test<-merge_data_batch_MYCN_normal[,-ind]
  set.seed(123)
  ind <- MYCN_conditions_new_amp$MYCN_conditions %>% createDataPartition(p = 0.7, list = FALSE)
  merge_data_batch_MYCN_amp_train<-merge_data_batch_MYCN_amp[,ind]
  merge_data_batch_MYCN_amp_test<-merge_data_batch_MYCN_amp[,-ind]
}
merge_data_batch_MYCN_normal_train<-merge_data_batch_MYCN_normal[,rownames(MYCN_conditions_new_normal_train)]
merge_data_batch_MYCN_normal_test<-merge_data_batch_MYCN_normal[,rownames(MYCN_conditions_new_normal_test)]

#### Consensus Cluster Plus ####
library("ConsensusClusterPlus")
dim(merge_data_batch_MYCN_normal_train)
d<-as.matrix(merge_data_batch_MYCN_normal_train)
11583*0.5 ### Train,Test gene number

mads=apply(d,1,mad) ### measured by median absolute deviation
d=d[rev(order(mads))[1:5792],] ### top 5000 variable genes # top 25%
d = sweep(d,1, apply(d,1,median,na.rm=T))
#### Figure 1b ####
results<- ConsensusClusterPlus(d,maxK=5,reps=1000,pItem=0.8,pFeature=1,title="Train_KM", ### Use This one
                               distance="euclidean",clusterAlg="km",seed=123,plot="pdf")

icl = calcICL(results,title="Train_KM",plot="pdf")
{Kvec = 2:
    x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
  PAC = rep(NA,length(Kvec))
  names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
  for(i in Kvec){
    M = results[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)}#end for i# The optimal K
  optK = Kvec[which.min(PAC)]}### help to determine the best K

a<-results[3][[1]]$consensusClass
a<-as.data.frame(a)
a$Sample<-rownames(a)
a_1<-a[a[,1]==1,] # Subgroup1
train_1<-merge_data_batch_MYCN_normal_train[,rownames(a_1)]
a_2<-a[a[,1]==2,] # Subgroup2
train_2<-merge_data_batch_MYCN_normal_train[,rownames(a_2)]
a_3<-a[a[,1]==3,] # Subgroup3
train_3<-merge_data_batch_MYCN_normal_train[,rownames(a_3)]
# d<-as.matrix(merge_data_batch_MYCN_normal_test) #same for test data



#### Figure 1c ####
pdf('PCA_removebatch_MYCN_Nomral_grey95.pdf',width = 8,height = 8)
fviz_pca_ind(pca_batch, col.ind=MYCN_conditions_new$Subgroup, mean.point=F, addEllipses = F, legend.title="Subgroup",
             geom.ind = "point")+scale_shape_manual(values=seq(0,13))+theme_classic()+
  scale_color_manual(values=c("#F8766D", "#7CAE00", "#00BFC4","#C77CFF","gray95"))
dev.off()


