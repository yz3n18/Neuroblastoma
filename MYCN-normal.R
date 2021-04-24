{discovery_datasets <- list( dataObj1,dataObj2,dataObj3,dataObj4,
                            dataObj5,dataObj6,dataObj7,dataObj8,dataObj9,
                            dataObj10,dataObj11,dataObj12,dataObj13,dataObj14,
                            dataObj15,dataObj16,dataObj17,dataObj18,dataObj19,
                            dataObj20)

names(discovery_datasets) = c( dataObj1$formattedName,
                               dataObj2$formattedNam,dataObj3$formattedName,
                               dataObj4$formattedName,dataObj5$formattedName,
                               dataObj6$formattedName,dataObj7$formattedName,
                               dataObj8$formattedName,dataObj9$formattedName,
                               dataObj10$formattedName,dataObj11$formattedName,
                               dataObj12$formattedName,dataObj13$formattedName,
                               dataObj14$formattedName,dataObj15$formattedName,
                               dataObj16$formattedName,dataObj17$formattedName,
                               dataObj18$formattedName,dataObj19$formattedName,dataObj20$formattedName)

validation_datasets <- list(dataObj21)
names(validation_datasets)<-c(dataObj21$formattedName)

exampleMetaObj=list() 
exampleMetaObj$originalData <- discovery_datasets
exampleMetaObj_validation<-list()
exampleMetaObj_validation$originalData <- validation_datasets
rm(merge_data_batch,merge_data_batch_MYCN_amp_test,merge_data_batch_MYCN_amp_train,merge_data_batch_MYCN_normal,merge_data_batch_MYCN_normal_test,
   merge_data_batch_MYCN_normal_train)
checkDataObject(exampleMetaObj, "Meta", "Pre-Analysis") ###Check your metaObject before MetaAnalysis using 
checkDataObject(exampleMetaObj_validation, "Meta", "Pre-Analysis")
exampleMetaObj <- runMetaAnalysis(exampleMetaObj, maxCores=4,runLeaveOneOutAnalysis = FALSE)
exampleMetaObj <- filterGenes(exampleMetaObj, isLeaveOneOut = F, FDRThresh = 0.05,effectSizeThresh = 1.3)
forwardRes <- forwardSearch( metaObject = exampleMetaObj, 
                             filterObject =exampleMetaObj$filterResults[[1]])

rocPlot(datasetObject = exampleMetaObj_validation$originalData$`E-MEXP_3517_GPL5175_42`,
        filterObject = forwardRes)} ### leave one out


a<-datasetObject$originalData$GSE3960$keys
a<-as.data.frame(a)

library(cluster)
library(factoextra)
sil <- silhouette(cutree(spe.ch.ward, k = k), spe.ch)
plot(results_test[4][[1]]$consensusTree)

for (i in 2:5){
  i<-3
  sil<-silhouette(Test_kmean[i][[1]]$consensusClass, dist(t(d)))
  fviz_silhouette(sil)+scale_colour_manual(values = c("#7CAE00", "#00BFC4", "#C77CFF"))
  ggsave(paste0('silhouette_Test_hc_pearson_wd_',i,'.pdf'),width = 6,height = 4,
         onefile=FALSE)

}

sil_pos<-sil[, , drop = FALSE]
rownames(sil_pos)<-colnames(merge_data_batch_MYCN_normal_test)
pos_sil_index <- which(sil_pos[, "sil_width"] > 0)
sil_pos<-sil_pos[pos_sil_index, , drop = FALSE]
colnames(train_2_sil)
write(colnames(target_RNA_3),'TARGET_RNA_group3.txt')

setwd('/Users/lefan/Rcode/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/')
sil_pos_train<-sil_pos
sil_pos_test<-sil_pos
a<-results_train[3][[1]]$consensusClass
a<-Train_kmean[3][[1]]$consensusClass

a<-as.data.frame(sil_pos_train)
a<-as.data.frame(a)
a$Sample<-rownames(a)
a_1<-a[a[,1]==1,]
train_1<-merge_data_batch_MYCN_amp_train[,rownames(a_1)]
a_2<-a[a[,1]==2,]
train_2<-merge_data_batch_MYCN_amp_train[,rownames(a_2)]
a_3<-a[a[,1]==3,]
train_3<-merge_data_batch_MYCN_amp_train[,rownames(a_3)]
a_4<-a[a[,1]==4,]
train_4<-merge_data_batch_MYCN_amp_train[,rownames(a_4)]
a_5<-a[a[,1]==5,]
train_5<-merge_data_batch_MYCN_amp_train[,rownames(a_5)]

a_1<-a[a[,1]==1,]
train_1<-merge_data_batch_MYCN_normal_train[,rownames(a_1)]
a_2<-a[a[,1]==2,]
train_2<-merge_data_batch_MYCN_normal_train[,rownames(a_2)]
a_3<-a[a[,1]==3,]
train_3<-merge_data_batch_MYCN_normal_train[,rownames(a_3)]
a_4<-a[a[,1]==4,]
train_4<-merge_data_batch_MYCN_normal_train[,rownames(a_4)]
a_5<-a[a[,1]==5,]
train_5<-merge_data_batch_MYCN_amp_train[,rownames(a_5)]



a<-results_test[3][[1]]$consensusClass
a<-Test_kmean[3][[1]]$consensusClass
a<-as.data.frame(sil_pos_test)
a<-as.data.frame(a)
a$Sample<-rownames(a)
a_1<-a[a[,1]==1,]
test_1<-merge_data_batch_MYCN_amp_test[,rownames(a_1)]
a_2<-a[a[,1]==2,]
test_2<-merge_data_batch_MYCN_amp_test[,rownames(a_2)]
a_3<-a[a[,1]==3,]
test_3<-merge_data_batch_MYCN_amp_test[,rownames(a_3)]
a_4<-a[a[,1]==4,]
test_4<-merge_data_batch_MYCN_amp_test[,rownames(a_4)]
a_5<-a[a[,1]==5,]
test_5<-merge_data_batch_MYCN_amp_test[,rownames(a_5)]


a_1<-a[a[,1]==1,]
test_1<-merge_data_batch_MYCN_normal_test[,rownames(a_1)]
a_2<-a[a[,1]==2,]
test_2<-merge_data_batch_MYCN_normal_test[,rownames(a_2)]
a_3<-a[a[,1]==3,]
test_3<-merge_data_batch_MYCN_normal_test[,rownames(a_3)]
a_4<-a[a[,1]==4,]
test_4<-merge_data_batch_MYCN_normal_test[,rownames(a_4)]
a_5<-a[a[,1]==5,]
train_5<-merge_data_batch_MYCN_amp_train[,rownames(a_5)]
#### data_distribution ####

table(sapply(colnames(test_3),function(x) strsplit(x,"_")[[1]][1]))

#### Survival ####
BiocManager::install('survminer')
library(survminer)

library(survival)
survival_data<-read.csv('/Users/yihuawang/Neuroblastoma/Hu/Survival_ALL.csv')
rownames(survival_data)<-survival_data$Sample

Clinical_test_normal<-na.omit(Clinical_test_normal)
survival_data<-survival_data[c(rownames(Clinical_test_normal)[Clinical_test_normal$High.risk=='No'],
                                      colnames(merge_data_batch_MYCN_amp_train),
                                      colnames(merge_data_batch_MYCN_amp_test)),]



survival_data<-na.omit(survival_data)
#survival_data$Group<-'MYCN-control'
#survival_data$Group<-'MYCN-AMP'
survival_data$Group<-NA

survival_data[colnames(train_1_sil),4]<-'Subgroup1' ## 6是根据Group所在的列
survival_data[colnames(train_2_sil),4]<-'Subgroup2' ## 6是根据Group所在的列
survival_data[colnames(train_3_sil),4]<-'Subgroup3' ## 6是根据Group所在的列
survival_data[colnames(merge_data_batch_MYCN_amp_train),4]<-'MYCN-AMP' ## 6是根据Group所在的列

survival_data[colnames(test_1_sil),4]<-'Subgroup1' ## 6是根据Group所在的
survival_data[colnames(test_2_sil),4]<-'Subgroup2' ## 6是根据Group所在的列
survival_data[colnames(test_3_sil),4]<-'Subgroup3' ## 6是根据Group所在的列
survival_data[colnames(merge_data_batch_MYCN_amp_test),4]<-'MYCN-AMP' ## 6是根据Group所在的列
survival_data<-na.omit(survival_data)
#colnames(survival_data)
#survival_data_AMP<-survival_data[survival_data[,4]!='MYCN-control',]
survival_data_AMP<-survival_data[survival_data[,4]!='MYCN-AMP',]
#survival_data_AMP<-survival_data
fit<-survfit(Surv(Day_Survival,Event) ~ Group, data=survival_data_AMP)


KMsurvival_plot<-ggsurvplot(fit,data=survival_data_AMP,pval = TRUE, #show p-value of log-rank test，显示log-rank分析得到的P值
                            #conf.int = TRUE, #添加置信区间
                            pval.size=10,
                            legend.labs =  c("Subgroup1",'Subgroup2','Subgroup3'),
                            
                            legend.title='', 
                            xlab = "Time in Day",   ###  customize X axis label.x的标time in years
                            #xlim=c(0,50),
                            break.x.by=1000, ###改变坐标轴间距
                            ylab=paste0('Overall survival Days'),
                            #ylab=paste0('Overall peak day'),
                            #ylab=paste0('SRAS−CoV−2 RNA +'),
                            surv.median.line = "hv", #添加中位生存时间的线
                            palette = c("#7CAE00", "#00BFC4", "#C77CFF"), ###  自定义颜色
                            #font.main = c(16, "bold", "darkblue"),
                            font.x = 40, # X 轴
                            font.y = 35, # c(14, "bold.italic", "darkred"), y 轴特征
                            font.tickslab = 30,# c(12, "plain", "darkgreen"), 坐标轴大小
                            #conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样子
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)的形式展示risk table
                            risk.table.y.text.col = T,###  colour risk table text annotations.
                            risk.table.y.text = F,###  show bars instead of names in text annotations in legend of risk table.不显示legend名字
                            tables.y.text=F,
                            tables.x.text =F,
                            #risk.table.title="My title", ## 改变 title名字
                            fontsize=5, ## 表格
                            ncensor.plot = F, #我TRUE就显示
                            #tables.theme=theme_cleantable(), # 取消table边框
                            ggtheme = theme_classic()#绘图主题
)
KMsurvival_plot$plot<-KMsurvival_plot$plot+
  theme(legend.text = element_text(size = 30),plot.margin = unit(c(2,2,0,2), "cm"))

KMsurvival_plot$table<-KMsurvival_plot$table+labs(x = NULL, y = NULL)+theme(axis.text=element_text(size=30),
                                                                            axis.title=element_text(size=30),
                                                                            
                                                                            legend.text = element_text(size = 25),
                                                                            plot.title = element_text(size=25),
                                                                            plot.margin = unit(c(0,2,2,2), "cm")) # 改变numer at risk 大小)

ggsave(file=paste0('All_silhouette_Km_3_survival_Lowrisk.pdf'),width = 15,height = 12,print(KMsurvival_plot),
       onefile=FALSE)
#### Clinical feature distribution ####
table(Clinical_test$Gender)

{Clinical_test_Amp<-Clinical_test[colnames(merge_data_batch_MYCN_amp),]
  Clinical_test_Amp<-na.omit(Clinical_test_Amp)
  table(Clinical_test_Amp$INSS.Stage)
  table(Clinical_test_Amp$Gender)
  table(Clinical_test_Amp$Age)
  table(Clinical_test_Amp$High.risk)
  ### Normal
  Clinical_test_normal<-Clinical_test[colnames(merge_data_batch_MYCN_normal),]
  Clinical_test_normal<-Clinical_test[colnames(merge_data_batch_MYCN_normal_train),]
  Clinical_test_normal<-Clinical_test[colnames(merge_data_batch_MYCN_normal_test),]
  Clinical_test_normal$Subgroup<-NA
  
  Clinical_test_normal[colnames(train_1),7]<-'1'
  Clinical_test_normal[colnames(train_2),7]<-'2'
  Clinical_test_normal[colnames(train_3),7]<-'3'
  
  Clinical_test_normal[colnames(test_1),7]<-'1'
  Clinical_test_normal[colnames(test_2),7]<-'2'
  Clinical_test_normal[colnames(test_3),7]<-'3'
  Clinical_test_normal<-na.omit(Clinical_test_normal)
  
  #INSS
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='1',]$INSS.Stage)
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='2',]$INSS.Stage)
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='3',]$INSS.Stage)
  
  #Gender
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='1',]$Gender)
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='2',]$Gender)
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='3',]$Gender)
  #Age
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='1',]$Age)
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='2',]$Age)
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='3',]$Age)
  #Risk
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='1',]$High.risk)
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='2',]$High.risk)
  table(Clinical_test_normal[Clinical_test_normal$Subgroup=='3',]$High.risk)
  
  
  survival_data_stage4<-survival_data[c(rownames(Clinical_test_normal)[Clinical_test_normal$High.risk=='No'],
                                        colnames(merge_data_batch_MYCN_amp_train),
                                        colnames(merge_data_batch_MYCN_amp_test)),]
  survival_data<-survival_data[c(rownames(Clinical_test_normal)[Clinical_test_normal$High.risk=='Yes']),]
  rownames(Clinical_test_Amp)[1:10]}

#### DEG ####
{train_1<-as.data.frame(train_1)
  dim(train_1)
  train_2<-as.data.frame(train_2)
  dim(train_2)
  train_3<-as.data.frame(train_3)
  dim(train_3)
  #train_4<-as.data.frame(train_4)
  #dim(train_4)
  
  merge_data_batch_MYCN_amp<-as.data.frame(merge_data_batch_MYCN_amp)
  dim(merge_data_batch_MYCN_amp)
  merge_data_batch_MYCN_amp_train_train1<-data.frame(merge_data_batch_MYCN_amp,
                                                        train_1,train_2,train_3)
  
  colnames(merge_data_batch_MYCN_amp_train_train1)<-c(paste0('AMP_',rep(1:714)),
                                                         paste0('Normal1_',rep(1:length(train_1[1,]))),
                                                         paste0('Normal2_',rep(1:length(train_2[1,]))),
                                                         paste0('Normal3_',rep(1:length(train_3[1,]))))
  rep(1,10)
  design <- model.matrix(~-1+factor(c(rep(1,714),rep(2,length(train_1[1,])),
                                      rep(3,length(train_2[2,])),rep(4,length(train_3[3,])))))
  colnames(design) <- c("AMP","Normal1","Normal2",'Normal3') # ??????design????????????rename the colnames of design
  
  contrastmatrix <- makeContrasts(Normal1-AMP,Normal2-AMP,Normal3-AMP,levels=design)##and make the contrasts
  
  
  fit <- lmFit(as.matrix(merge_data_batch_MYCN_amp_train_train1), design) # Run limma to get Differentially expressed genes
  fit2 <- contrasts.fit(fit, contrastmatrix)  # Run limma to get Differentially expressed genes
  ?lmFit
  #the borrowed variance approach described in class
  fit2 <- eBayes(fit2) # ??????????????????FDR????????????
  myresults_train_1<-topTable(fit2,coef=1, # 1 for FIH, 2 for VHL, 3 for VHL_FIH
                              adjust="fdr",
                              number=nrow(merge_data_batch_MYCN_amp_train_train1))
  myresults_train1_UP<-subset(myresults_train_1, adj.P.Val < 0.001)
  myresults_train1_UP_500<-myresults_train1_UP[order(myresults_train1_UP$logFC,decreasing = T), ][1:500,]
  myresults_train1_DOWN_500<-myresults_train1_UP[order(myresults_train1_UP$logFC,decreasing = F), ][1:500,]
  myresults_train_2<-topTable(fit2,coef=2, # 1 for FIH, 2 for VHL, 3 for VHL_FIH
                              adjust="fdr",
                              number=nrow(merge_data_batch_MYCN_amp_train_train1))
  
  myresults_train2_UP<-subset(myresults_train_2, adj.P.Val < 0.001)
  myresults_train2_UP_500<-myresults_train2_UP[order(myresults_train2_UP$logFC,decreasing = T), ][1:500,]
  myresults_train2_DOWN_500<-myresults_train2_UP[order(myresults_train2_UP$logFC,decreasing = F), ][1:500,]
  
  
  myresults_train_3<-topTable(fit2,coef=3, # 1 for FIH, 2 for VHL, 3 for VHL_FIH
                              adjust="fdr",
                              number=nrow(merge_data_batch_MYCN_amp_train_train1))
  myresults_train3_UP<-subset(myresults_train_3, adj.P.Val < 0.001)
  myresults_train3_UP_500<-myresults_train3_UP[order(myresults_train3_UP$logFC,decreasing = T), ][1:500,]
  myresults_train3_DOWN_500<-myresults_train3_UP[order(myresults_train3_UP$logFC,decreasing = F), ][1:500,]
  
  myresults_train_4<-topTable(fit2,coef=4, # 1 for FIH, 2 for VHL, 3 for VHL_FIH
                              adjust="fdr",
                              number=nrow(merge_data_batch_MYCN_amp_train_train1))
  myresults_train4_UP<-subset(myresults_train_4, adj.P.Val < 0.001)
  myresults_train4_UP_500<-myresults_train4_UP[order(myresults_train4_UP$logFC,decreasing = T), ][1:500,]
  myresults_train4_DOWN_500<-myresults_train4_UP[order(myresults_train4_UP$logFC,decreasing = F), ][1:500,]} #### train

 d_train_FUDAN<-as.data.frame(d_train_FUDAN)

write.table(myresults_train_2,'MYCN-Nomral-Train2+AMP.txt',sep = '\t',row.names =T )
write.table(merge_data_batch_MYCN_amp_train_train1,'MYCN-Nomral-Train+AMP.txt',sep = '\t',row.names =T )
?write.table
{test_1<-as.data.frame(test_1)
  dim(test_1)
  test_2<-as.data.frame(test_2)
  dim(test_2)
  test_3<-as.data.frame(test_3)
  dim(test_3)
  test_4<-as.data.frame(test_4)
  dim(test_4)
  
  merge_data_batch_MYCN_amp<-as.data.frame(merge_data_batch_MYCN_amp)
  dim(merge_data_batch_MYCN_amp)
  merge_data_batch_MYCN_amp_test_test1<-data.frame(merge_data_batch_MYCN_amp,
                                                   test_1,test_2,test_3)
  
  colnames(merge_data_batch_MYCN_amp_test_test1)<-c(paste0('AMP_',rep(1:714)),
                                                    paste0('Normal1_',rep(1:length(test_1[1,]))),
                                                    paste0('Normal2_',rep(1:length(test_2[1,]))),
                                                    paste0('Normal3_',rep(1:length(test_3[1,]))))
  rep(1,10)
  design <- model.matrix(~-1+factor(c(rep(1,714),rep(2,length(test_1[1,])),
                                      rep(3,length(test_2[1,])),rep(4,length(test_3[1,])))))
  colnames(design) <- c("AMP","Normal1","Normal2",'Normal3') # ??????design????????????rename the colnames of design

  contrastmatrix <- makeContrasts(Normal1-AMP,Normal2-AMP,Normal3-AMP,levels=design)##and make the contrasts
  
  
  fit <- lmFit(as.matrix(merge_data_batch_MYCN_amp_test_test1), design) # Run limma to get Differentially expressed genes
  fit2 <- contrasts.fit(fit, contrastmatrix)  # Run limma to get Differentially expressed genes
  #the borrowed variance approach described in class
  fit2 <- eBayes(fit2) # ??????????????????FDR????????????
  myresults_test_1<-topTable(fit2,coef=1, # 1 for FIH, 2 for VHL, 3 for VHL_FIH
                             adjust="fdr",
                             number=nrow(merge_data_batch_MYCN_amp_test_test1))
  myresults_test1_UP<-subset(myresults_test_1, adj.P.Val < 0.001)
  myresults_test1_UP_500<-myresults_test1_UP[order(myresults_test1_UP$logFC,decreasing = T), ][1:500,]
  myresults_test1_DOWN_500<-myresults_test1_UP[order(myresults_test1_UP$logFC,decreasing = F), ][1:500,]
  myresults_test_2<-topTable(fit2,coef=2, # 1 for FIH, 2 for VHL, 3 for VHL_FIH
                             adjust="fdr",
                             number=nrow(merge_data_batch_MYCN_amp_test_test1))
  
  myresults_test2_UP<-subset(myresults_test_2, adj.P.Val < 0.001)
  myresults_test2_UP_500<-myresults_test2_UP[order(myresults_test2_UP$logFC,decreasing = T), ][1:500,]
  myresults_test2_DOWN_500<-myresults_test2_UP[order(myresults_test2_UP$logFC,decreasing = F), ][1:500,]
  
  
  myresults_test_3<-topTable(fit2,coef=3, # 1 for FIH, 2 for VHL, 3 for VHL_FIH
                             adjust="fdr",
                             number=nrow(merge_data_batch_MYCN_amp_test_test1))
  myresults_test3_UP<-subset(myresults_test_3, adj.P.Val < 0.001)
  myresults_test3_UP_500<-myresults_test3_UP[order(myresults_test3_UP$logFC,decreasing = T), ][1:500,]
  myresults_test3_DOWN_500<-myresults_test3_UP[order(myresults_test3_UP$logFC,decreasing = F), ][1:500,]
  
  myresults_test_4<-topTable(fit2,coef=4, # 1 for FIH, 2 for VHL, 3 for VHL_FIH
                             adjust="fdr",
                             number=nrow(merge_data_batch_MYCN_amp_test_test1))
  myresults_test4_UP<-subset(myresults_test_4, adj.P.Val < 0.001)
  myresults_test4_UP_500<-myresults_test4_UP[order(myresults_test4_UP$logFC,decreasing = T), ][1:500,]
  myresults_test4_DOWN_500<-myresults_test4_UP[order(myresults_test4_UP$logFC,decreasing = F), ][1:500,]} #### test
rm(myresults_1,myresults_1_DOWN_500,myresults_1_UP,myresults_1_UP_500,
   myresults_test1_UP_500,myresults_test1_DOWN_500,myresults_train1_UP_500,myresults_train1_DOWN_500,
   myresults_2,myresults_2_DOWN_500,myresults_2_UP,myresults_2_UP_500,
   myresults_test2_UP_500,myresults_test2_DOWN_500,myresults_train2_UP_500,myresults_train2_DOWN_500,
   myresults_3,myresults_3_DOWN_500,myresults_3_UP,myresults_3_UP_500,
   myresults_test3_UP_500,myresults_test3_DOWN_500,myresults_train3_UP_500,myresults_train3_DOWN_500,
   myresults_4,myresults_4_DOWN_500,myresults_4_UP,myresults_4_UP_500,
   myresults_test4_UP_500,myresults_test4_DOWN_500,myresults_train4_UP_500,myresults_train4_DOWN_500)
rm(myresults_test_1,myresults_test1_UP,myresults_train_1,myresults_train1_UP,
   myresults_test_2,myresults_test2_UP,myresults_train_2,myresults_train2_UP,
   myresults_test_3,myresults_test3_UP,myresults_train_3,myresults_train3_UP,
   myresults_test_4,myresults_test4_UP,myresults_train_4,myresults_train4_UP)

length(intersect(rownames(myresults_test1_UP_500),rownames(myresults_train1_UP_500)))
length(intersect(rownames(myresults_test2_UP_500),rownames(myresults_train1_UP_500)))
length(intersect(rownames(myresults_test3_UP_500),rownames(myresults_train1_UP_500)))
#length(intersect(rownames(myresults_test4_UP_500),rownames(myresults_train4_UP_500)))

length(intersect(rownames(myresults_test1_DOWN_500),rownames(myresults_train1_DOWN_500)))
length(intersect(rownames(myresults_test2_DOWN_500),rownames(myresults_train1_DOWN_500)))
length(intersect(rownames(myresults_test3_DOWN_500),rownames(myresults_train1_DOWN_500)))
#length(intersect(rownames(myresults_test4_DOWN_500),rownames(myresults_train4_DOWN_500)))

length(intersect(rownames(myresults_test1_UP_500),rownames(myresults_train2_UP_500)))
length(intersect(rownames(myresults_test2_UP_500),rownames(myresults_train2_UP_500)))
length(intersect(rownames(myresults_test3_UP_500),rownames(myresults_train2_UP_500)))
#length(intersect(rownames(myresults_test4_UP_500),rownames(myresults_train4_UP_500)))

length(intersect(rownames(myresults_test1_DOWN_500),rownames(myresults_train2_DOWN_500)))
length(intersect(rownames(myresults_test2_DOWN_500),rownames(myresults_train2_DOWN_500)))
length(intersect(rownames(myresults_test3_DOWN_500),rownames(myresults_train2_DOWN_500)))
#length(intersect(rownames(myresults_test4_DOWN_500),rownames(myresults_train4_DOWN_500)))

length(intersect(rownames(myresults_test1_UP_500),rownames(myresults_train3_UP_500)))
length(intersect(rownames(myresults_test2_UP_500),rownames(myresults_train3_UP_500)))
length(intersect(rownames(myresults_test3_UP_500),rownames(myresults_train3_UP_500)))
#length(intersect(rownames(myresults_test4_UP_500),rownames(myresults_train4_UP_500)))

length(intersect(rownames(myresults_test1_DOWN_500),rownames(myresults_train3_DOWN_500)))
length(intersect(rownames(myresults_test2_DOWN_500),rownames(myresults_train3_DOWN_500)))
length(intersect(rownames(myresults_test3_DOWN_500),rownames(myresults_train3_DOWN_500)))
#length(intersect(rownames(myresults_test4_DOWN_500),rownames(myresults_train4_DOWN_500)))

#### 临床瀑布图 ####

Clinical_test_Amp<-Clinical_test[colnames(merge_data_batch_MYCN_amp),]
Clinical_test_Amp$Subgroup<-'1'
Clinical_test_Amp[colnames(train_3),7]<-'3'
Clinical_test_Amp[colnames(test_3),7]<-'1'
Clinical_test_Amp$Gene1<-1
Clinical_test_Amp$Gene2<-2

Clinical_test_normal_t<-t(Clinical_test_Amp)

Clinical_test_normal<-Clinical_test[colnames(merge_data_batch_MYCN_normal),]
Clinical_test_normal<-Clinical_test[colnames(merge_data_batch_MYCN_normal_train),]

Clinical_test_normal$Subgroup<-NA
Clinical_test_normal[colnames(train_1),7]<-'1'
Clinical_test_normal[colnames(train_2),7]<-'2'
Clinical_test_normal[colnames(train_3),7]<-'3'
Clinical_test_normal[colnames(train_4),7]<-'4'

Clinical_test_normal[colnames(test_1),7]<-'2'
Clinical_test_normal[colnames(test_2),7]<-'4'
Clinical_test_normal[colnames(test_3),7]<-'3'
Clinical_test_normal[colnames(test_4),7]<-'1'

Clinical_test_normal$Gene1<-1
Clinical_test_normal$Gene2<-2
Clinical_test_normal_t<-t(Clinical_test_normal)


Clinical_test_normal_t<-Clinical_test_normal_t[,order(Clinical_test_normal_t[7,])]
dim(Clinical_test_Amp_t)
dim(Clinical_test_normal_t)

Clinical_test_normal_t_89<-as.data.frame(Clinical_test_normal_t[8:9,])
Clinical_test_normal_t_89[1:2,]<-sapply(Clinical_test_normal_t_89[1:2,], as.numeric)

pdf('pdf.pdf',width = 100,height = 3)
ha = HeatmapAnnotation(Age118months = Clinical_test_normal_t[2,],
                       INSSstage =Clinical_test_normal_t[3,],
                       MYCNaaamplification=Clinical_test_normal_t[4,],
                       Male=Clinical_test_normal_t[5,],
                       HighRisk=Clinical_test_normal_t[6,],
                       Subgroup=Clinical_test_normal_t[7,],
                       border = F,
                       gp = gpar(col = "white"),
                       #annotation_label = c("Age >18 months", "INSS stage", "MYCN amplification",
                       #                    'Male','High risk','Subgroup'),
                       col = list(Age118months = c("No" = "grey", "Yes" = "Black", "Not defined" = "azure2"),
                                  INSSstage = c("1" = "thistle1", "2" = "hotpink", "3" = "deeppink",
                                                '4'='hotpink4','4S'='goldenrod1',"Not defined" = "azure2"),
                                  MYCNaaamplification = c("No" = "grey", "Yes" = "Black", "Not defined" = "azure2"),
                                  Male = c("No" = "grey", "Yes" = "Black", "Not defined" = "azure2"),
                                  HighRisk = c("No" = "grey", "Yes" = "Black", "Not defined" = "azure2"),
                                  Subgroup = c("1" = "cadetblue1", "2" = "dodgerblue3", "3" = "plum3",'4'='darkorchid4')
                       ),
                       annotation_legend_param = list(
                         Age118months = list(
                           title = "Age >18 months",
                           direction = "horizontal",nrow = 1),
                         INSSstage = list(
                           title = "INSS stage",
                           direction = "horizontal",nrow = 1),
                         MYCNaaamplification = list(
                           title = "MYCN amplification",
                           direction = "horizontal",nrow = 1),
                         Male = list(
                           title = "Male",
                           direction = "horizontal",nrow = 1),
                         HighRisk = list(
                           title = "High risk",
                           direction = "horizontal",nrow = 1),
                         Subgroup = list(
                           title = "Subgroup",
                           direction = "horizontal",nrow = 1)
                         
                       ),
                       annotation_name_side = "left",
                       annotation_height=unit(10, "cm"))
draw(Heatmap(as.matrix(Clinical_test_normal_t_89[1,]),
             column_labels=rep(' ',length(Clinical_test_normal_t_89[1,])),
             bottom_annotation = ha,
             heatmap_legend_param = list(direction = "horizontal")),
     annotation_legend_side = "bottom")

dev.off()
table(Clinical_test_normal[Clinical_test_normal$Subgroup=='1',]$High.risk)
table(Clinical_test_normal[Clinical_test_normal$Subgroup=='2',]$High.risk)
table(Clinical_test_normal[Clinical_test_normal$Subgroup=='3',]$High.risk)
table(Clinical_test_normal[Clinical_test_normal$Subgroup=='4',]$High.risk)


X<-matrix(c(275,194,126,98,67,
            54,54,89,383,63,
            185,174,94,362,59,110,49,32,35,141),nrow=5)
X<-matrix(c(235,166,91,68,57,
            22,30,59,259,47,
            124,108,56,255,41,
            72,36,20,21,78),nrow=5)
X<-matrix(c(40,28,35,30,10,
            32,24,30,124,16,
            61,66,38,107,18,
            38,13,12,14,63),nrow=5)
X<-matrix(c(135,199,
            174,103,
            189,190,
            51,123),nrow=2)
X<-matrix(c(103,168,
            116,59,
            124,132,
            32,76),nrow=2)
X<-matrix(c(32,31,
            58,44,
            65,58,
            19,47),nrow=2)
X<-matrix(c(257,482,
            364,255,
            442,404,
            80,279),nrow=2)
X<-matrix(c(198,401,
            251,149,
            312,259,
            43,180),nrow=2)
X<-matrix(c(59,81,
            113,106,
            130,145,
            37,99),nrow=2)
X<-matrix(c(16,147,
            77,57,
            68,106,
            3,77),nrow=2)
X<-matrix(c(8,122,
            53,35,
            48,63,
            1,46),nrow=2)
X<-matrix(c(8,25,
            24,22,
            20,43,
            2,31),nrow=2)
X
chisq.test(X)


#### Subtype siganture ####
train_1<-as.data.frame(train_1)
dim(train_1)
train_2<-as.data.frame(train_2)
dim(train_2)
train_3<-as.data.frame(train_3)
dim(train_3)
train_4<-as.data.frame(train_4)
dim(train_4)
1034+370
merge_data_batch_MYCN_amp_train<-as.data.frame(merge_data_batch_MYCN_amp_train)
dim(merge_data_batch_MYCN_amp_train)
rm(merge_data_batch_MYCN_amp_train_train1)
merge_data_batch_MYCN_amp_train_train1<-data.frame(train_2,train_1,train_3)
dim(merge_data_batch_MYCN_amp_train_train1)
colnames(merge_data_batch_MYCN_amp_train_train1)<-c(paste0('Subgroup2_',rep(1:316)),
                                                    paste0('nonSubgroup2_',rep(1:1404)))
write.table(merge_data_batch_MYCN_amp_train_train1,'Subgroup2vsnon.txt',sep = '\t') 
write.table(t(c(rep('Subgroup2',316),rep('nonSubgroup2',1404))),'Subgroup1_name.txt',sep = ' ')


#### GSEA火山图 ####
library(xlsx)
GSEA_2vs1<-read.csv('GSEA_1vs2.csv',header = T,sep = ',',stringsAsFactors = F)
GSEA_2vs3<-read.csv('GSEA_1vs3.csv',header = T,sep = ',',stringsAsFactors = F)
GSEA_2vs13<-merge(GSEA_2vs1,GSEA_2vs3,by='Pathway')
write.csv(GSEA_2vs13,'GSEA_1vs23.csv')
GSEA_2vs13$finalgroup<-'NA'
for (i in 1:length(GSEA_2vs13[,1])){
  if (GSEA_2vs13[i,5] == GSEA_2vs13[i,9]){
    GSEA_2vs13[i,10]<-GSEA_2vs13[i,5]
  } else {
    GSEA_2vs13[i,10]<-'Confused'
  }
}
require(dplyr)
library(ggrepel)
p<-GSEA_2vs13 %>%
  mutate(bin = Sub1vs3NES < 1) %>%
  ggplot(aes(x=Sub1vs2NES,y=Sub1vs3NES,color =Dataset.x))+
  scale_color_manual(values=c("gray28",'brown3', "cadetblue3",'goldenrod1','thistle3'))+## change the color manually
  geom_point(alpha=0.7,size = 3)+# change transparency
  #ylim(0,15)+# change range of x/y axis
  #scale_x_continuous(limits=c(1.2,2))+
  #scale_y_continuous(limits=c(-2,2.5))+
  labs(x='\nNES-Subgroup [1/2]',y='NES-Subgroup [1/3]\n')+### add title and x/y label
  theme_classic(base_family = "Avenir")+##change type of theme
  facet_grid(bin ~ Dataset.x, scale='free_y')+
  geom_smooth(method='lm',formula= y~x,se = FALSE)+
  #geom_text_repel(data=subset(GSEA_2vs13,FDR.x< 0.05))+
                  
  #geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+## add line on vertical
  #geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+#+ add line on horizontal
  theme(legend.position="right",## legend is annotation of group
        panel.grid=element_blank(),
        strip.text.y = element_blank(),
        strip.text.x =element_text(size=14, color="black",
                     face="bold"),
        legend.title = element_blank(),### delete the title of legend
        legend.text= element_text(face="bold", color="black",family = "Times", size=15),##set font
        axis.text.x = element_text(face="bold", color="black", family = "Times",size=15),
        axis.text.y = element_text(face="bold", color="black",family = "Times", size=15),
        axis.title.x = element_text(face="bold", color="black",family = "Times", size=20),
        axis.title.y = element_text(face="bold",color="black",family = "Times", size=20),
        plot.margin = unit(c(1,1,1,1), "cm"))
ggsave("volcano_plot of ras_tgfb.png", p, dpi = 300, width = 20, height = 8)
ggsave(p,filename = "volcano_plot of ras_tgfb.pdf", width = 16, height = 8)
ggsave("volcano_plot of ras_tgfb.pdf", p,  width = 10, height = 8)
?strip.text.y
#colnames(merge_data_batch_MYCN_amp_train_train1)[1:10]

merge_data_batch_MYCN_amp_train_train1<-data.frame(train_1,train_2,train_3)
matrix_name<-matrix(data = NA,nrow = 1720,ncol = 2)
matrix_name[1:1034,1]<-paste0(rep('Subgroup1_',1034),1:1034)
matrix_name[1035:1350,1]<-paste0(rep('Subgroup2_',316),1:316)
matrix_name[1351:1720,1]<-paste0(rep('Subgroup3_',370),1:370)
write.csv(matrix_name,'matrix_name.csv')
rm(matrix_name)
write.table(merge_data_batch_MYCN_amp_train_train1,'Train_group.txt',sep = '\t')
dim(train_1)# 1034
dim(train_2)# 316
dim(train_3) # 370
dim(merge_data_batch_MYCN_normal_test_estimate) # 378
dim(test_2) # 108
dim(test_3) # 164
379+108+164
matrix_name<-matrix(data = NA,nrow = 650,ncol = 2)
matrix_name[1:378,1]<-paste0(rep('Subgroup1_',378),1:378)
matrix_name[379:486,1]<-paste0(rep('Subgroup2_',108),1:108)
matrix_name[487:650,1]<-paste0(rep('Subgroup3_',164),1:164)
write.csv(matrix_name,'matrix_name.csv')
write.table(merge_data_batch_MYCN_normal_test_estimate,'Test_group.txt',sep = '\t')

gene<-'TERT'
median(train_1[gene,])
median(train_2[gene,])
median(train_3[gene,])
median(as.matrix(merge_data_batch_MYCN_amp)[gene,])

median(test_1['TERT',])
median(test_2['AURKA',])
median(test_3['AURKA',])
316+370
1034+316
median(as.matrix(merge_data_batch_MYCN_amp['ODC1',]))

colnames(merge_data_batch_MYCN_amp_train_train1)<-c(paste0('Subgroup1_',rep(1:1034)),
                                                    paste0('Subgroup2_',rep(1:316)),
                                                    paste0('Subgroup2_',rep(1:370))
                                                    )
colnames(merge_data_batch_MYCN_amp_train_train1)<-c(paste0('Subgroup2_',rep(1:370)),
                                                    paste0('Subgroup1_',rep(1:1350))
)

dim(merge_data_batch_MYCN_amp_train_train1)

design <- model.matrix(~-1+factor(c(rep(1,1034),rep(2,316),rep(3,370))))
design <- model.matrix(~-1+factor(c(rep(1,370),rep(2,1350))))
colnames(design) <- c("Subgroup1",'Subgroup2','Subgroup3') # ??????design????????????rename the colnames of design
colnames(design) <- c("Subgroup2",'Subgroup1')
contrastmatrix <- makeContrasts(Subgroup2-Subgroup1,Subgroup2-Subgroup3,levels=design)##and make the contrasts
contrastmatrix <- makeContrasts(Subgroup2-Subgroup1,levels=design)

fit <- lmFit(as.matrix(merge_data_batch_MYCN_amp_train_train1), design) # Run limma to get Differentially expressed genes
fit2 <- contrasts.fit(fit, contrastmatrix)  # Run limma to get Differentially expressed genes
fit2 <- eBayes(fit2) # ??????????????????FDR????????????
limma_train_subgroup2_results<-topTable(fit2, coef=1,adjust="fdr",
                                        number=nrow(merge_data_batch_MYCN_amp_train_train1))
limma_train_subgroup2_results$Gene<-rownames(limma_train_subgroup2_results)
limma_train_subgroup1_results<-merge(limma_train_subgroup1_results,limma_train_subgroup2_results,by='Gene')

write.table(limma_train_subgroup1_results,'IPA_Subgroup1.txt',sep = '\t')
limma_train_subgroup1_results_Up<-subset(limma_train_subgroup1_results,c(logFC > 0.5849625 & adj.P.Val < 0.001))
limma_train_subgroup1_results_down<-subset(limma_train_subgroup1_results,c(logFC < -0.5849625 & adj.P.Val < 0.001))
#write.table(limma_train_subgroup1_results,'limma_train_subgroup1vs2_results_FDR.txt',row.names=T,sep = '\t')
write.table(limma_train_subgroup1_results_Up,'limma_train_subgroup2vs1_results_Up_FDR_FC1.5.txt',row.names=T,sep = '\t')
write.table(limma_train_subgroup1_results_down,'limma_train_subgroup2vs1_results_Down_FDR_FC1.5.txt',row.names=T,sep = '\t')
log2(1.5)
2^1
#### Venn plot for overlapping DEG ####
limma_train_subgroup1_results_Up<-read.table('limma_train_subgroup2vs1_results_Up_FDR.txt',
                                           header = T,sep='\t')
rownames(limma_train_subgroup1_results_Up)<-limma_train_subgroup1_results_Up$Gene
limma_train_subgroup1_results_down<-read.table('limma_train_subgroup2vs1_results_Down_FDR.txt',
                                             header = T,sep='\t')
rownames(limma_train_subgroup1_results_down)<-limma_train_subgroup1_results_down$Gene
limma_train_subgroup1_results_Up2<-read.table('limma_train_subgroup2vs3_results_Up_FDR.txt',
                                             header = T,sep='\t')
rownames(limma_train_subgroup1_results_Up2)<-limma_train_subgroup1_results_Up2$Gene
limma_train_subgroup1_results_dow2<-read.table('limma_train_subgroup2vs3_results_Down_FDR.txt',
                                               header = T,sep='\t')
rownames(limma_train_subgroup1_results_dow2)<-limma_train_subgroup1_results_dow2$Gene
venn.diagram(list(Sub2vs1_Up=as.factor(limma_train_subgroup1_results_Up$Gene),
                  Sub2vs1_Down=as.factor(limma_train_subgroup1_results_down$Gene),
                  Sub2vs3_Up=as.factor(limma_train_subgroup1_results_Up2$Gene),
                  Sub2vs3_Down=as.factor(limma_train_subgroup1_results_dow2$Gene)),
             resolution = 300, imagetype = "tiff", 
             col='transparent',
             fill=c("lightblue3", "mediumpurple2", "red3",'goldenrod1'),
             main="Differentially expressed proteins",
             main.cex = 3, main.fontface = 2,alpha = 0.50, cex=2.5, cat.cex=2.5,
             margin=0.2, #改变图片大小
             cat.dist = c(0.24, 0.24, 0.1,0.1),##注释距离圆的距离
             # cex for each area label
             #cat.cex for category name
             #cat.col = c("darkblue", "darkgreen", "orange"),
             #print.mode = 'percent',
             #category.names = c('RAS(75)',expression(paste('TGF-',beta,"(35)")),expression(paste('RAS+TGF-',beta,"(47)"))),
             filename = "Subgroup2_DEG_FDR_FC1.5.tiff")
subgroup1up<-intersect(limma_train_subgroup1_results_Up$Gene,
                                     limma_train_subgroup1_results_Up2$Gene)
limma_train_subgroup1_results_Up<-limma_train_subgroup1_results_Up[subgroup1up,]
limma_train_subgroup1_results_Up2<-limma_train_subgroup1_results_Up2[subgroup1up,]

subgroup1up<-intersect(limma_train_subgroup1_results_down$Gene,
                       limma_train_subgroup1_results_dow2$Gene)
limma_train_subgroup1_results_Up<-limma_train_subgroup1_results_down[subgroup1up,]
limma_train_subgroup1_results_Up2<-limma_train_subgroup1_results_dow2[subgroup1up,]

limma_train_subgroup1_results_Up<-merge(limma_train_subgroup1_results_Up,
                                         limma_train_subgroup1_results_Up2,by='Gene')

limma_train_subgroup1_results_Up<-limma_train_subgroup1_results_Up[,c(1,2,6,8,12)]
write.csv(limma_train_subgroup1_results_Up,'Subgroup2_DEG_FDR.csv')
getwd()


#### fgsea ####
{library(fgsea)
?fgsea
data(examplePathways)
data(exampleRanks)
fgseaRes <- fgsea(pathways = examplePathways, 
                  stats = exampleRanks,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)
head(fgseaRes[order(pval), ])
sum(fgseaRes[, padj < 0.01])
plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
               exampleRanks) + labs(title="Programmed Cell Death"))
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam = 0.5)

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
                                      examplePathways, exampleRanks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(examplePathways[mainPathways], exampleRanks, fgseaRes, 
              gseaParam = 0.5)
library(data.table)
fwrite(fgseaRes, file="fgseaRes.txt", sep="\t", sep2=c("", " ", ""))}
#### GSVA ####
library(GSVA)

library(fgsea)
gmt_file <- "/Users/lefan/Rcode/Neuroblastoma/Hu/GEO+Array/Metabolism.gmt"
pathways <- gmtPathways(gmt_file)
gs <- lapply(pathways, unique)
# 接下来去掉那pathways里出现过的genes
count <- table(unlist(gs))
keep <- names(which(table(unlist(gs)) < 2))
gs <- lapply(gs, function(x) intersect(keep, x))
# 过滤之后，很多pathway一个gene都不剩了，去掉这些
gs <- gs[lapply(gs, length) > 0]
gsva_R <- gsva(as.matrix(merge_data_batch_MYCN_amp_train_train1), pathways,method='ssgsea')
dim(gsva_R)
colnames(gsva_R)<-c(paste0('Normal1_',rep(1:731)),
                    paste0('Normal2_',rep(1:476)),
                    paste0('Normal3_',rep(1:677)),
                    paste0('Normal4_',rep(1:276)))

  
design <- model.matrix(~-1+factor(c(rep(1,731),rep(2,476),
                                    rep(3,677),rep(4,276))))

colnames(design) <- c("Subgroup1",'Subgroup2',
                      'Subgroup3','Subgroup4') # ??????design????????????rename the colnames of design

contrastmatrix <- makeContrasts(Subgroup1-Subgroup2,Subgroup1-Subgroup3,
                                Subgroup1-Subgroup4,Subgroup2-Subgroup3,
                                Subgroup2-Subgroup4,Subgroup3-Subgroup4,levels=design)##and make the contrasts

gene<-"MYC"
mean(c(as.matrix(train_2[gene,]),as.matrix(train_1[gene,]),as.matrix(train_3[gene,])))
mean(as.matrix(train_2[gene,]))
mean(as.matrix(train_2[gene,]))

mean(as.matrix(merge_data_batch_MYCN_amp_train[gene,]))

colnames(train_1)
fit <- lmFit(as.matrix(gsva_R), design) # Run limma to get Differentially expressed genes
fit2 <- contrasts.fit(fit, contrastmatrix)  # Run limma to get Differentially expressed genes
fit2 <- eBayes(fit2) # ??????????????????FDR????????????
gsva_R_result<-topTable(fit2, adjust="fdr",
                                        number=nrow(gsva_R))
gsva_R_result<-matrix(nrow=50,ncol=4)
rownames(gsva_R_result)<-rownames(gsva_R)
colnames(gsva_R_result)<-c("Subgroup1",'Subgroup2',
                           'Subgroup3','Subgroup4')

for (i in 1:50){
  gsva_R_result[i,1]<-mean(gsva_R[i,1:731])
  gsva_R_result[i,2]<-mean(gsva_R[i,732:1207])
  gsva_R_result[i,3]<-mean(gsva_R[i,1208:1884])
  gsva_R_result[i,4]<-mean(gsva_R[i,1885:2160])
  
}
dim(gsva_R)
write.table(gsva_R_result,'gsva_Hallmark_result.txt',sep = '\t')
### one-way anova
gsva_R_t<-matrix(nrow = 50,ncol = 2)

for (i in 1:50){
  a<-shapiro.test(gsva_R[i,1:731])
  b<-shapiro.test(gsva_R[i,732:1207])
  c<-shapiro.test(gsva_R[i,1208:1884])
  d<-shapiro.test(gsva_R[i,1885:2160])
  if (min(a$p.value,b$p.value,c$p.value,d$p.value) < 0.05) {
    d<-kruskal.test(list(gsva_R[i,1:731],gsva_R[i,732:1207],
                         gsva_R[i,1208:1884],gsva_R[i,1885:2160]))
    gsva_R_t[i,1]<-rownames(gsva_R)[i]
    gsva_R_t[i,2]<-d$p.value
  }
  
}
class(gsva_R[1,2160])
?kruskal.test
 
#print('u-test')
#print(c$p.value)

result[1,i]<-'Kruskal-Wallis'
result[2,i]<-c$p.value



#### Tumor purity ####
library(utils)
rforge <- "http://r-forge.r-project.org"
install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(stringr)
?estimateScore
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

?write.table

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
?estimateScore

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
#Purity_All[intersect(rownames(Purity_All),colnames(train_4)),6]<-'Subgroup_4'
Purity_All[intersect(rownames(Purity_All),colnames(test_1)),5]<-'Subgroup_1'
Purity_All[intersect(rownames(Purity_All),colnames(test_2)),5]<-'Subgroup_2'
Purity_All[intersect(rownames(Purity_All),colnames(test_3)),5]<-'Subgroup_3'
#Purity_All[intersect(rownames(Purity_All),colnames(test_4)),6]<-'Subgroup_1'

dim(merge_data_batch_MYCN_normal_train)

dim(train_2)
setdiff(colnames(train_2),rownames(Purity_Train[Purity_Train$Group=='Subgroup_2',]))
setdiff(rownames(Purity_Train[Purity_Train$Group=='Subgroup_2',]),colnames(train_2))
Purity_Train<-Purity_All[c(colnames(merge_data_batch_MYCN_normal_train),
                           colnames(merge_data_batch_MYCN_amp_train)),]
merge_data_batch_MYCN_normal_train_estimate<-data.frame(train_1,train_2,train_3)

Purity_Train<-Purity_All[c(intersect(rownames(Purity_All),colnames(merge_data_batch_MYCN_normal_train_estimate)),
                           intersect(rownames(Purity_All),colnames(merge_data_batch_MYCN_amp))),]
write.csv(merge_data_batch_MYCN_normal_train_estimate[1:5,],'tets.csv')
merge_data_batch_MYCN_normal_train_estimate[1,'EMTAB8248_US22502540_252038210022_S01_GE1-v5_10_Apr08_1_2']
colnames(train_2)
Purity_Train<-na.omit(Purity_Train)
Purity_Train$Group<-'MYCN_AMP'
Purity_Train[intersect(rownames(Purity_Train),colnames(train_1)),6]<-'Subgroup_1'
Purity_Train[intersect(rownames(Purity_Train),colnames(train_2)),6]<-'Subgroup_2'
Purity_Train[intersect(rownames(Purity_Train),colnames(train_3)),6]<-'Subgroup_3'
Purity_Train[intersect(rownames(Purity_Train),colnames(train_4)),6]<-'Subgroup_4'
Purity_Train[colnames(train_1),6]<-'Subgroup_1'
Purity_Train<-na.omit(Purity_Train)
Purity_Train[colnames(train_2),6]<-'Subgroup_2'
Purity_Train<-na.omit(Purity_Train)
Purity_Train[colnames(train_3),6]<-'Subgroup_3'
Purity_Train<-na.omit(Purity_Train)
#Purity_Train[colnames(train_4),5]<-'Subgroup_4'
dim(Purity_Train)

Purity_Test<-Purity_All[c(colnames(merge_data_batch_MYCN_normal_test),
                          colnames(merge_data_batch_MYCN_amp_test)),]
merge_data_batch_MYCN_normal_test_estimate<-data.frame(test_1,test_2,test_3)
Purity_Test<-Purity_All[c(intersect(rownames(Purity_All),colnames(merge_data_batch_MYCN_normal_test_estimate)),
                           intersect(rownames(Purity_All),colnames(merge_data_batch_MYCN_amp))),]


Purity_Test<-na.omit(Purity_Test)

Purity_Test$Group<-'MYCN_AMP'
Purity_Test[colnames(test_1),6]<-'Subgroup_1'
Purity_Test<-na.omit(Purity_Test)
Purity_Test[colnames(test_2),6]<-'Subgroup_2'
Purity_Test<-na.omit(Purity_Test)
Purity_Test[colnames(test_3),6]<-'Subgroup_3'
Purity_Test<-na.omit(Purity_Test)
#Purity_Test[colnames(test_4),6]<-'Subgroup1'
#Purity_Test<-na.omit(Purity_Test)

rownames(Purity_Train)<-Purity_Train$Sample
fonts()
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
229+387
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

log(0)
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


#### Complex Heatmap & WGCNA ####
### K-mean 3
{a<-Train_kmean[3][[1]]$consensusClass
a<-as.data.frame(a)
a$Sample<-rownames(a)


a_1<-a[a[,1]==1,]
train_1<-merge_data_batch_MYCN_normal_train[,rownames(a_1)]
length(intersect(rownames(a_1),colnames(merge_data_batch_MYCN_normal_train)))

a_2<-a[a[,1]==2,]
train_2<-merge_data_batch_MYCN_normal_train[,rownames(a_2)]
a_3<-a[a[,1]==3,]
train_3<-merge_data_batch_MYCN_normal_train[,rownames(a_3)]
a_4<-a[a[,1]==4,]
train_4<-merge_data_batch_MYCN_normal_train[,rownames(a_4)]
a_5<-a[a[,1]==5,]
train_5<-merge_data_batch_MYCN_amp_train[,rownames(a_5)]

dim(merge_data_batch_MYCN_normal_train)

a<-results_test[3][[1]]$consensusClass
a<-as.data.frame(a)
a$Sample<-rownames(a)


a_1<-a[a[,1]==1,]
test_1<-merge_data_batch_MYCN_normal_test[,rownames(a_1)]
a_2<-a[a[,1]==2,]
test_2<-merge_data_batch_MYCN_normal_test[,rownames(a_2)]
a_3<-a[a[,1]==3,]
test_3<-merge_data_batch_MYCN_normal_test[,rownames(a_3)]
a_4<-a[a[,1]==4,]
test_4<-merge_data_batch_MYCN_normal_test[,rownames(a_4)]
a_5<-a[a[,1]==5,]
train_5<-merge_data_batch_MYCN_amp_train[,rownames(a_5)]}

merge_data_batch_MYCN_amp_train_train1
d<-as.matrix(merge_data_batch_MYCN_normal_train)
11583*0.5
mads=apply(d,1,mad) ### measured by median absolute deviation
d=d[rev(order(mads))[1:5792],] ### top 5000 variable genes # top 25%


# pheatmap(pheatmap(y[cg,],show_rownames = F,show_colnames = F),scale = "row")
merge_data_batch_MYCN_amp_train_train1<-data.frame(merge_data_batch_MYCN_amp,train_1,train_2,train_3)
merge_data_batch_MYCN_amp_train_train1<-data.frame(train_1,train_2,train_3)
merge_data_batch_MYCN_amp_train_train1_wgcna<-t(merge_data_batch_MYCN_amp_train_train1[rownames(d),])
merge_data_batch_MYCN_amp_amp_wgcna<-t(data.frame(merge_data_batch_MYCN_amp)[rownames(d),])

dim(merge_data_batch_MYCN_amp_train_train1_wgcna)
dim(merge_data_batch_MYCN_amp_amp_wgcna)
n=t(scale(merge_data_batch_MYCN_amp_train_train1_wgcna)) # 'scale'可以对log-ratio
n=t(scale(merge_data_batch_MYCN_amp_amp_wgcna)) # 'scale'可以对log-ratio数值进行归一化

n[n>3]=3
n[n< -3]= -3
n<-n[,1:714]
dim(n)

n[1:4,1:4]
options(stringsAsFactors = F)
ac<-data.frame(c(rep('AMP',714),rep('Subgroup1',1034),rep('Subgroup2',316), # AMP+Train
                 rep('Subgroup3',370)))
ac<-data.frame(c(rep('Subgroup1',1034),rep('Subgroup2',316),# Train only
                 rep('Subgroup3',370)))
ac<-data.frame(c(rep('AMP',714))) # AMP only
colnames(ac)<-'Molecular subtype'
rownames(ac)=colnames(n) #把ac的行名给到n的列名，即对每一个探针标记上分组信息'noTNBC'还是'TNBC'）
#ac$`Molecular subtype` <- factor(ac$`Molecular subtype`, levels = c("AMP", "Subgroup1",'Subgroup2','Subgroup3'))
ac$Blank<-'1'
ac<-as.matrix(t(ac))
colorbar<-colorRampPalette(c('royalblue3','gray92','firebrick2'))(n=1000)
library(ComplexHeatmap)
set.seed(123)
pdf("heatmap.pdf", width = 10, height = 10)
ha = HeatmapAnnotation(`Molecular subtype` =ac[1,],
                       col = list(`Molecular subtype` = c('AMP' = "#F8766D", 'Subgroup1' = "#7CAE00",
                                                          'Subgroup2'='#00BFC4','Subgroup3'='#C77CFF')),
                       annotation_legend_param = list(
                         `Molecular subtype` = list(title = "Molecular subtype")))
draw(Heatmap(n, name = "Value", column_km = 4,column_km_repeats = 1000,
             clustering_distance_columns = "euclidean",
             show_row_names = FALSE,show_column_names  = FALSE,
             top_annotation = ha,col=colorbar,show_heatmap_legend = T,
             show_parent_dend_line = FALSE),padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()
set.seed(123)
pdf("heatmap_train_normal_order_coldend.pdf", width = 10, height = 10)
ha = HeatmapAnnotation(`Molecular subtype` =ac[1,],
                       col = list(`Molecular subtype` = c('Subgroup1' = "#7CAE00",
                                                          'Subgroup2'='#00BFC4','Subgroup3'='#C77CFF')),
                       annotation_legend_param = list(
                         `Molecular subtype` = list(title = "Molecular subtype")))
ht<-draw(Heatmap(n, name = "Value", column_km = 3,column_km_repeats = 1000,
                 clustering_distance_columns = "euclidean",
                 show_row_names = FALSE,show_column_names  = FALSE,
                 top_annotation = ha,col=colorbar,show_heatmap_legend = T,
                 column_dend_reorder = TRUE,
                 show_parent_dend_line = FALSE),padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()
row_order_ht<-row_order(ht)
row_dend(ht)
merge_data_batch_MYCN_amp_amp_wgcna<-merge_data_batch_MYCN_amp_amp_wgcna[,row_order_ht]
n<-n[row_order_ht,]
pdf("heatmap_train_amp_order4.pdf", width = 5, height = 9.7)
ha = HeatmapAnnotation(`Molecular subtype` =ac[1,],
                       col = list(`Molecular subtype` = c('AMP' = "#F8766D")),
                       annotation_legend_param = list(
                         `Molecular subtype` = list(title = "Molecular subtype")))
draw(Heatmap(n, name = "Value", column_km = 1,column_km_repeats = 1,
             clustering_distance_columns = "euclidean",
             show_row_names = FALSE,show_column_names  = FALSE,
             top_annotation = ha,col=colorbar,show_heatmap_legend = T,
             cluster_rows = F,
             show_parent_dend_line = FALSE),padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()
rownames_n<-as.data.frame(rownames(n))

n %>%
  rownames()%>%
  as.data.frame()%>%
  write.csv('heatmap_rownames.csv')


?Heatmap
merge_data_batch_MYCN_amp_train_train1_wgcna<-t(merge_data_batch_MYCN_amp_train_train1[rownames(d),]) ### top 5000 variable genes # top 25%
dim(merge_data_batch_MYCN_amp_train_train1_wgcna)



{
  library(WGCNA)
  library(doParallel) ### 多线程运行
  registerDoParallel(cores=4)
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  sft = pickSoftThreshold(merge_data_batch_MYCN_amp_train_train1_wgcna, 
                          powerVector = powers, verbose = 5)
  ?pickSoftThreshold
  #设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
  # Plot the results:
  ##sizeGrWindow(9, 5)
  
  # Scale-free topology fit index as a function of the soft-thresholding power
  pdf(file = paste("SFTPlot.pdf",sep = "."),width = 10,height = 7)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  # this line corresponds to using an R^2 cut-off of h
  #abline(h=0.90,col="red")
  abline(h=0.85,col="blue")
  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  dev.off()
  sft$powerEstimate
  #### 提取模块中基因信息 ####
  # Select module
  module = "brown";
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
  inModule = (moduleColors==module);
  modProbes = probes[inModule];
  
  #### WGCNA step by step ####
  datExpr<-merge_data_batch_MYCN_amp_train_train1_wgcna
  softPower = 5
  adjacency = adjacency(datExpr, power = softPower,type='signed hybrid')
  ?adjacency
  ## Topological Overlap Matrix (TOM)
  TOM = TOMsimilarity(adjacency,
                      TOMType = "signed")
  dissTOM = 1-TOM
  ## Clustering using TOM
  geneTree = hclust(as.dist(dissTOM), method = "average")#
  # Plot the resulting clustering tree (dendrogram)
  pdf('WGCNA_geneTree_Signhybrid.pdf',width = 100,height = 10)
  plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
       labels = FALSE, hang = 0.04)
  dev.off()
  # We like large modules, so we set the minimum module size relatively high:
  minModuleSize = 100
  # Module identification using dynamic tree cut:
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minModuleSize);
  table(dynamicMods)
  # Convert numeric lables into colors
  dynamicColors = labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  sizeGrWindow(8,6)
  pdf('WGCNA_geneTree_50.pdf',width = 100,height = 10)
  plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module colors")
  dev.off()
  # Calculate eigengenes
  MEList = moduleEigengenes(datExpr, colors = dynamicColors)
  MEs = MEList$eigengenes
  MEDiss = 1-cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average")
  MEDissThres = 0.25
  plotEigengeneNetworks(MEs, 
                        "Eigengene adjacency heatmap", 
                        marHeatmap = c(3,4,2,2), 
                        plotDendrograms = FALSE, 
                        xLabelsAngle = 90) 
  
  
  # Calculate dissimilarity of module eigengenes
  # Plot the result
  pdf('WGCNA_MergedgeneTree_50.pdf',width = 100,height = 10)
  plot(METree, main = "Clustering of module eigengenes",
       xlab = "", sub = "")
  abline(h=MEDissThres, col = "red")
  dev.off()
  
  # Call an automatic merging function
  merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs
  table(mergedColors)
  pdf('WGCNA_MergedgeneTree_new_50.pdf',width = 10,height = 10)
  plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                      c("Dynamic Tree Cut", "Merged dynamic"),
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  
  ####临床-WGCNA ####
  {
    design=model.matrix(~0+ as.factor(datTraits$Subtype))
    colnames(design)=levels(as.factor(datTraits$Subtype))
    moduleColors <- mergedColors
    # Recalculate MEs with color labels
    MEs0 = moduleEigengenes(merge_data_batch_MYCN_amp_train_train1_wgcna, 
                            moduleColors)$eigengenes
    MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs)
    moduleTraitCor = cor(MEs, design , use = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
    
    sizeGrWindow(10,6)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    pdf("WGCNA_phenotype_200.pdf", width = 10, height = 10)
    par(mar = c(6, 8.5, 3, 3));
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(design),
                   yLabels = names(MEs),
                   ySymbols = names(MEs),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,
                   zlim = c(-1,1),
                   main = paste("Module-trait relationships"))
    dev.off()
  }
  # 感兴趣性状的模块的具体基因分析
  modNames = substring(names(MEs), 3) # MEturquoise从第三个开始取成，turquoise
  geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))# 基因和颜色模块的correlation
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), 
                                            nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="");
  
  #再计算性状与基因的相关性矩阵
  ## 只有连续型性状才能只有计算
  ## 这里把 CB 表型这个变量用0,1进行数值化。
  CB = as.data.frame(design[,3]) 
  names(CB) = "CB"
  geneTraitSignificance = as.data.frame(cor(datExpr, CB, use = "p"));
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 
                                            nSamples));
  names(geneTraitSignificance) = paste("GS.", names(CB), sep="");
  names(GSPvalue) = paste("p.GS.", names(CB), sep="")
  
  #两个相关性矩阵联合起来,指定感兴趣模块进行分析
  
  par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
  par(mar = c(4,5,3,1))
  pdf('Module membership vs. gene significance.pdf',
      width = 20,height = 10)
  
  for (i in modNames){
    module = i
    column = match(i, modNames);
    moduleGenes = moduleColors==module;# 提取module中的gene
    verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                       abs(geneTraitSignificance[moduleGenes, 1]),
                       xlab = paste("Module Membership in", module, "module"),
                       ylab = "Gene significance for Subgroup3",
                       main = paste("Module membership vs gene significance\n"),
                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, 
                       col = module,
                       pch=19)
  }
  dev.off()
  
  # 提取指定模块的基因名
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基
  # 创建empty data frame
  module_gene<-as.data.frame(matrix(data=NA,nrow=max(table(moduleColors)),
                                    ncol=length(table(moduleColors))))
  colnames(module_gene)<-modNames
  # Select module
  for (i in 1:length(module_gene[1,])){
    module<-modNames[i]
    inModule = (moduleColors==module)
    j<-as.numeric(table(inModule)[2])
    print(j)
    module_gene[1:j,i]<-probes[inModule]
  }
  write.csv(module_gene,'module_gene_200.csv',row.names = F,quote = F)
  
  ####导出到cytoscape ####
  # Recalculate topological overlap
  TOM = TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate,
                              networkType='signed hybrid',
                              TOMType='signed') 
  # Select module
  module = "turquoise";
  # Select module probes
  probes = colnames(datExpr) ## 我们例子里面的probe就是基因名
  inModule = (moduleColors==module);
  modProbes = probes[inModule]; 
  ## 也是提取指定模块的基因名
  # Select the corresponding Topological Overlap
  modTOM = TOM[inModule, inModule];
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(
    modTOM,
    edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
    nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
    weighted = TRUE,
    threshold = 0.3,
    nodeNames = modProbes, 
    nodeAttr = moduleColors[inModule]
  )
  ?exportNetworkToCytoscape
  
  ####提取hub genes ####
  #计算连
  connet=abs(cor(datExpr,use="p"))^6
  Alldegrees1=intramodularConnectivity(connet, moduleColors)
  head(Alldegrees1)
  # 模块内的连通性与gene显著性的关系
  # (2) Relationship between gene significance and intramodular connectivity
  which.module="turquoise"
  EB= as.data.frame(design[,2]); # change specific 
  names(EB) = "EB"
  GS1 = as.numeric(cor(EB,datExpr, use="p"))
  GeneSignificance=abs(GS1)
  colorlevels=unique(moduleColors)
  sizeGrWindow(9,6)
  par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
  par(mar = c(4,5,3,1))
  for (i in c(1:length(colorlevels)))
  {
    whichmodule=colorlevels[[i]];
    restrict1 = (moduleColors==whichmodule);
    verboseScatterplot(Alldegrees1$kWithin[restrict1],
                       GeneSignificance[restrict1], col=moduleColors[restrict1],
                       main=whichmodule,pch=19,
                       xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
  }
  #计算模块内所有基因的连通性, 筛选hub gene
  ###(3) Generalizing intramodular connectivity for all genes on the array
  datKME=signedKME(datExpr, MEs, outputColumnName="MM.")
  # Display the first few rows of the data frame
  head(datKME)
  ##Finding genes with high gene significance and high intramodular connectivity in
  # interesting modules
  # abs(GS1)> .9 可以根据实际情况调整参数
  # abs(datKME$MM.black)>.8 至少大于 >0.8
  FilterGenes= abs(GS1)> .5 & abs(datKME$MM.turquoise)>.8
  table(FilterGenes)
  probes[FilterGenes]
}

#### gephi ####
merge_data_batch_MYCN_amp_train_train1_checkpoint<-read.table('/Users/lefan/Rcode/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/IPA_Subgroup1.txt',
                                                              header = T,sep = '\t',row.names = 1)
merge_data_batch_MYCN_amp_train_train1_checkpoint['MMP11',]
WGCNA_modulegene<-read.csv('/Users/lefan/Rcode/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/WGCNA_100.csv',
                           header = T,sep = ',') ## for gephi network construction
merge_data_batch_MYCN_amp_train_train1_checkpoint$Gene<-rownames(merge_data_batch_MYCN_amp_train_train1_checkpoint)
gephi_PPI<-merge(merge_data_batch_MYCN_amp_train_train1_checkpoint,
                 WGCNA_modulegene,by='Gene')
WGCNA_node<-read.csv('Test_node.csv',header = T,sep = ',')
WGCNA_node<-merge(WGCNA_node,gephi_PPI,by='Gene')
write.csv(WGCNA_node,'WGCNA_node.csv')
##### module genes with fold change
merge_data_batch_MYCN_amp_train_train1_checkpoint<-read.table('/Users/lefan/Rcode/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/IPA_Subgroup1.txt',
                                                              header = T,sep = '\t',row.names = 1)
WGCNA_modulegene<-read.csv('/Users/lefan/Rcode/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/WGCNA/module_gene_100.csv',
                           header = T,sep = ',') ## for immune checkpoint
rownames(merge_data_batch_MYCN_amp_train_train1_checkpoint)[1]
merge_data_batch_MYCN_amp_train_train1_checkpoint<-merge_data_batch_MYCN_amp_train_train1_checkpoint[WGCNA_modulegene$brown,]
dim(merge_data_batch_MYCN_amp_train_train1_checkpoint)
merge_data_batch_MYCN_amp_train_train1_checkpoint<-na.omit(merge_data_batch_MYCN_amp_train_train1_checkpoint)
write.csv(merge_data_batch_MYCN_amp_train_train1_checkpoint,'WGCNA_Subgroup3_brown_100.csv')


rm(GSPvalue,column,module,moduleGenes,geneTraitSignificance,MMPvalue,
   geneModuleMembership,modNames,moduleLabels,net,TOM,modTOM,dissTOM)
rm(FilterGenes,dynamicColors,dynamicMods,inModule,iSample,datExpr,
   datKME,EB,MEs,MEs0,METree,module_gene,moduleTraitPvalue)
rm(GS1,GeneSignificance,ModuleSignificance,which.trait,
   moduleTraitCor,y,ADJ1,Alldegrees1  )

#### Immune Signature+ cibersort+xcell+Timer #### 未成功

library(ggplot2)
Cibersort<-read.table('./Immune environment/CIBERSORT.Output_Job13.txt',header = T,sep = '\t')
Cibersort<-subset(Cibersort,P.value < 0.05)
Cibersort_1<-as.data.frame(matrix(data=NA,nrow = 22*1587,ncol = 4))
colnames(Cibersort_1)<-c("Sample", "Group", "Fraction", "Cell")
j<-0
for (i in 3:24){
  j<-j+1
  Cibersort_1[c(1+1587*(j-1)):(1587*j),1:3]<-Cibersort[,c(1,2,i)]
  Cibersort_1[c(1+1587*(j-1)):(1587*j),4]<-c(rep(colnames(Cibersort)[i],1587))
  
}

my_comparisons_train <- list( c("Subgroup1", "MYCN-AMP"), 
                              c("Subgroup2", "MYCN-AMP"),
                              c("Subgroup3", "MYCN-AMP"),
                              c('Subgroup1','Subgroup2'),
                              c('Subgroup1','Subgroup3'),
                              c('Subgroup2','Subgroup3'))
ggplot(Cibersort_1, aes(x=Cell, y=Fraction, fill=Group)) + 
  geom_boxplot()+theme_classic()+stat_compare_means(comparisons = my_comparisons_train,size=10)+
  theme(axis.text.x = element_text( face="bold",color = 'black',angle=45,hjust = 1))
ggsave('Cibersort.pdf',width = 40,height = 10)
1:2434
2435:4868
?matrix
dim(Cibersort)
library(ggplot2)
library(dplyr)
library(forcats)
library(hrbrthemes)
library(viridis)
Cibersort_1 %>%
  mutate(Cell = fct_reorder(Cell, Fraction)) %>%
  mutate(Cell = factor(Cell)) %>%
  ggplot(aes(fill=Group, y=Fraction, x=Cell)) + 
  geom_violin(position="dodge", alpha=0.5, outlier.colour="transparent") +
  scale_fill_viridis(discrete=T, name="") +
  theme_ipsum()  +
  xlab("") +
  ylab("Tip (%)") +
  ylim(0,1)


ggplot(Cibersort_1,aes(fill=  Group, y = Fraction,x = Cell))+geom_violin()+
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
###Xcell ####
# use website to get result
### the number of sample overpass Timer's limitation
#### Immune checkpoint ####
merge_data_batch_MYCN_amp_train_train1<-data.frame(train_1,train_2,train_3)
Immune_checkpoint<-read.csv('/Users/lefan/Rcode/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/Immune environment/Immune_checkpoint.csv',
                            header = T,sep = ',')
merge_data_batch_MYCN_amp_train_train1_checkpoint<-read.table('/Users/lefan/Rcode/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/IPA_Subgroup1.txt',
                                                              header = T,sep = '\t',row.names = 1)

colnames(merge_data_batch_MYCN_amp_train_train1_checkpoint)
merge_data_batch_MYCN_amp_train_train1_checkpoint<-merge_data_batch_MYCN_amp_train_train1_checkpoint[,c(1,5,6,10,11,15)]
merge_data_batch_MYCN_amp_train_train1_checkpoint<-merge_data_batch_MYCN_amp_train_train1_checkpoint[Immune_checkpoint[,2],]
merge_data_batch_MYCN_amp_train_train1_checkpoint$Group<-Immune_checkpoint[,3]
merge_data_batch_MYCN_amp_train_train1_checkpoint<-na.omit(merge_data_batch_MYCN_amp_train_train1_checkpoint)
merge_data_batch_MYCN_amp_train_train1_checkpoint_pvalue<-merge_data_batch_MYCN_amp_train_train1_checkpoint[,c(2,4,6)]
merge_data_batch_MYCN_amp_train_train1_checkpoint_pvalue$`P-value`<-apply(merge_data_batch_MYCN_amp_train_train1_checkpoint_pvalue,
                                                                          1,min)

pch = rep("na", 51)
is_sig = merge_data_batch_MYCN_amp_train_train1_checkpoint_pvalue$`P-value` < 0.05
pch[is_sig] = '*'
is_sig = merge_data_batch_MYCN_amp_train_train1_checkpoint_pvalue$`P-value` < 0.01
pch[is_sig] = '**'
is_sig = merge_data_batch_MYCN_amp_train_train1_checkpoint_pvalue$`P-value` < 0.001
pch[is_sig] = '***'
is_sig = merge_data_batch_MYCN_amp_train_train1_checkpoint_pvalue$`P-value` < 0.0001
pch[is_sig] = '****'

min(merge_data_batch_MYCN_amp_train_train1_checkpoint_pvalue)


write.csv(merge_data_batch_MYCN_amp_train_train1_checkpoint,'checkpoint.csv')
merge_data_batch_MYCN_amp_train_train1_checkpoint<-t(merge_data_batch_MYCN_amp_train_train1_checkpoint)
n_immunecheckpoint<-merge_data_batch_MYCN_amp_train_train1_checkpoint[1:3,]
n_immunecheckpoint<-apply(n_immunecheckpoint, 2,as.numeric)


rownames(n_immunecheckpoint)<-c('Subgroup1',
                                'Subgroup2',
                                'Subgroup3')
ac_immunecheckpoint<-as.data.frame(merge_data_batch_MYCN_amp_train_train1_checkpoint[4,])
colnames(ac_immunecheckpoint)<-'Immune checkpoint'

rownames(ac_immunecheckpoint)[2]<-'Immune checkpoint'
ac_immunecheckpoint$`Immune checkpoint` <- factor(ac_immunecheckpoint$`Immune checkpoint`, 
                                                  levels = c('None','Inhibitory','Stimulatory'))
ac_immunecheckpoint$Blank<-'1'
ac_immunecheckpoint<-as.matrix(t(ac_immunecheckpoint))
colorcheckpoint<-colorRampPalette(c('royalblue4','lightseagreen','gray92','tan2','firebrick3'))(n=1000)

pdf("heatmap_Immunecheckpoint_ordernew.pdf", width = 17, height = 3.4)
ha = HeatmapAnnotation(`P-value`=anno_text(pch),
  `Immune checkpoint` =ac_immunecheckpoint[1,],
                       col = list(`Immune checkpoint` = c( 'Inhibitory' = "sienna4",
                                               'None'='gray17',
                                               'Stimulatory' = "pink")),
                       annotation_legend_param = list(
                         `Immune checkpoint`= list(title = "Immune checkpoint")),
                       annotation_name_gp= gpar(fontsize = 15)
                       )
draw(Heatmap(n_immunecheckpoint,
             clustering_distance_columns = "euclidean",
             show_row_names = T,show_column_names  = T,
             column_names_rot = 75,
             row_names_gp = gpar(fontsize = 18),
             top_annotation = ha,col=colorcheckpoint,show_heatmap_legend = T,
             heatmap_legend_param = list(
               title = "logFC", at = c(-2,-1, 0,1, 2)),
             show_parent_dend_line = FALSE),padding = unit(c(1, 1, 1, 1), "cm"))
dev.off()
library(ComplexHeatmap)
### 热图，row是log1/2&3，column是gene


#### MYCN expression ####
my_comparisons_train <- list(c('MYCN-AMP','Subgroup2'),
                              c('Subgroup1','Subgroup2'),
                              c('Subgroup1','Subgroup3'),
                              c('Subgroup2','Subgroup3'))
train_1_sil['AURKA',1]
dim(merge_data_batch_MYCN_amp_train_train1)

merge_data_batch_MYCN_amp_train_train1<-data.frame(merge_data_batch_MYCN_amp,
                                                   train_1_sil,train_2_sil,train_3_sil)
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['AURKA',]
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['MDM2',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['ODC1',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['CCND1',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['AKT1',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['PIK3CA',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['CDK4',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['CDK6',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['PHOX2B',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['TWIST1',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['ASCL1',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['LMO1',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['MYC',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['TP53',] # 无mdm

merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['JAZF1',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['CENPA',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['HAMP',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['MYCN',] # 无mdm
merge_data_batch_MYCN_amp_train_train1<-merge_data_batch_MYCN_amp_train_train1['MYC',] # 无mdm


merge_data_batch_MYCN_amp_train_train1<-t(merge_data_batch_MYCN_amp_train_train1)
merge_data_batch_MYCN_amp_train_train1<-as.data.frame(merge_data_batch_MYCN_amp_train_train1)
merge_data_batch_MYCN_amp_train_train1$Group<-'NA'
dim(merge_data_batch_MYCN_amp_train_train1)

merge_data_batch_MYCN_amp_train_train1[1:714,2]<-'MYCN-AMP'
merge_data_batch_MYCN_amp_train_train1[715:1748,2]<-'Subgroup1'
merge_data_batch_MYCN_amp_train_train1[1749:2064,2]<-'Subgroup2'
merge_data_batch_MYCN_amp_train_train1[2065:2434,2]<-'Subgroup3'
merge_data_batch_MYCN_amp_train_train1<-na.omit(merge_data_batch_MYCN_amp_train_train1)
## palette = c("#7CAE00", "#00BFC4", "#C77CFF") 三组配色，加在add.params后面


ggviolin(merge_data_batch_MYCN_amp_train_train1, x = "Group", y = 'MYC',
         fill = "Group",
         add = "boxplot", add.params = list(fill = "white"),palette = c("#D22C6C", "#6688AB", "#97C17E",'#BA86B5'))+
  stat_compare_means(comparisons = my_comparisons_train,size=10)+
  ylab("MYC expression\n") +
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
ggsave('Subgroup_AMP_MYCN_expression.pdf',width = 13,height = 12) # 3组 宽12， 4组 宽 16
show_col(hue_pal()(4))
#### MYC-score ####
MYC_hallmark<-read.csv('MYC_hallmark.csv',header = T,sep = ',')
#MYC_hallmark<-list(MYC_hallmark$MYC1)
gsva_MYC_matrix<-as.data.frame(matrix(data=NA,
                                      nrow=3084,ncol=7))
colnames(gsva_MYC_matrix)<-c('AURKA','MYC1-score','MYC2-score','MYC-score',
                             'Group','Group-MYCN','Group-subgroup')

gsva_MYC_matrix[1:714,1]<-t(merge_data_batch_MYCN_amp['AURKA',])
gsva_MYC_matrix[1:714,2]<-t(gsva(as.matrix(merge_data_batch_MYCN_amp),list(MYC_hallmark$MYC1) , mx.diff=1,kcdf="Gaussian",
                                 method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[1:714,3]<-t(gsva(as.matrix(merge_data_batch_MYCN_amp),list(MYC_hallmark$MYC2) , mx.diff=1,kcdf="Gaussian",
                                 method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[1:714,4]<-t(gsva(as.matrix(merge_data_batch_MYCN_amp),list(MYC_hallmark$MYC) , mx.diff=1,kcdf="Gaussian",
                                 method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[1:714,5]<-'Neuroblastoma'
gsva_MYC_matrix[1:714,6]<-'MYCN-AMP'
gsva_MYC_matrix[1:714,7]<-'MYCN-AMP'

dim(train_1)
train_1_mean<-apply(train_3, 1,mean)
write.table(train_1_mean,'train_3_mean.txt',sep = '\t')

?apply
714+1034

gsva_MYC_matrix[715:1748,1]<-(train_1['AURKA',])
gsva_MYC_matrix[715:1748,2]<-t(gsva(as.matrix(train_1),list(MYC_hallmark$MYC1) , mx.diff=1,kcdf="Gaussian",
                                 method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[715:1748,3]<-t(gsva(as.matrix(train_1),list(MYC_hallmark$MYC2) , mx.diff=1,kcdf="Gaussian",
                                 method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[715:1748,4]<-t(gsva(as.matrix(train_1),list(MYC_hallmark$MYC) , mx.diff=1,kcdf="Gaussian",
                                 method= "ssgsea",ssgsea.norm=F))
dim(test_1)
1748+378


gsva_MYC_matrix[1749:2126,1]<-t(test_1['AURKA',])
gsva_MYC_matrix[1749:2126,2]<-t(gsva(as.matrix(test_1),list(MYC_hallmark$MYC1) , mx.diff=1,kcdf="Gaussian",
                                    method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[1749:2126,3]<-t(gsva(as.matrix(test_1),list(MYC_hallmark$MYC2) , mx.diff=1,kcdf="Gaussian",
                                    method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[1749:2126,4]<-t(gsva(as.matrix(test_1),list(MYC_hallmark$MYC) , mx.diff=1,kcdf="Gaussian",
                                    method= "ssgsea",ssgsea.norm=F))


gsva_MYC_matrix[715:2126,5]<-'Neuroblastoma'
gsva_MYC_matrix[715:2126,6]<-'MYCN-Normal'
gsva_MYC_matrix[715:2126,7]<-'Subgroup1'

dim(train_2)
2126+316
gsva_MYC_matrix[2127:2442,1]<-(train_2['AURKA',])
gsva_MYC_matrix[2127:2442,2]<-t(gsva(as.matrix(train_2),list(MYC_hallmark$MYC1) , mx.diff=1,kcdf="Gaussian",
                                    method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2127:2442,3]<-t(gsva(as.matrix(train_2),list(MYC_hallmark$MYC2) , mx.diff=1,kcdf="Gaussian",
                                    method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2127:2442,4]<-t(gsva(as.matrix(train_2),list(MYC_hallmark$MYC) , mx.diff=1,kcdf="Gaussian",
                                    method= "ssgsea",ssgsea.norm=F))
dim(test_2)
2442+108
gsva_MYC_matrix[2443:2550,1]<-t(test_2['AURKA',])
gsva_MYC_matrix[2443:2550,2]<-t(gsva(as.matrix(test_2),list(MYC_hallmark$MYC1) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2443:2550,3]<-t(gsva(as.matrix(test_2),list(MYC_hallmark$MYC2) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2443:2550,4]<-t(gsva(as.matrix(test_2),list(MYC_hallmark$MYC) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))


gsva_MYC_matrix[2127:2550,5]<-'Neuroblastoma'
gsva_MYC_matrix[2127:2550,6]<-'MYCN-Normal'
gsva_MYC_matrix[2127:2550,7]<-'Subgroup2'

dim(train_3)
2550+370
gsva_MYC_matrix[2551:2920,1]<-(train_3['AURKA',])
gsva_MYC_matrix[2551:2920,2]<-t(gsva(as.matrix(train_3),list(MYC_hallmark$MYC1) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2551:2920,3]<-t(gsva(as.matrix(train_3),list(MYC_hallmark$MYC2) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2551:2920,4]<-t(gsva(as.matrix(train_3),list(MYC_hallmark$MYC) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
dim(test_3)
2920+164
gsva_MYC_matrix[2921:3084,1]<-t(test_3['AURKA',])
gsva_MYC_matrix[2921:3084,2]<-t(gsva(as.matrix(test_3),list(MYC_hallmark$MYC1) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2921:3084,3]<-t(gsva(as.matrix(test_3),list(MYC_hallmark$MYC2) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2921:3084,4]<-t(gsva(as.matrix(test_3),list(MYC_hallmark$MYC) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))


gsva_MYC_matrix[2551:3084,5]<-'Neuroblastoma'
gsva_MYC_matrix[2551:3084,6]<-'MYCN-Normal'
gsva_MYC_matrix[2551:3084,7]<-'Subgroup3'

gsva_MYC_matrix[1:714,8]<-colnames(merge_data_batch_MYCN_amp)
gsva_MYC_matrix[715:1748,8]<-colnames(train_1)
gsva_MYC_matrix[1749:2126,8]<-colnames(test_1)
gsva_MYC_matrix[2127:2442,8]<-colnames(train_2)
gsva_MYC_matrix[2443:2550,8]<-colnames(test_2)
gsva_MYC_matrix[2551:2920,8]<-colnames(train_3)
gsva_MYC_matrix[2921:3084,8]<-colnames(test_3)

gsva_MYC_matrix[1:714,9]<-'MYCN-AMP'
gsva_MYC_matrix[715:1748,9]<-'Train1'
gsva_MYC_matrix[1749:2126,9]<-'Test1'
gsva_MYC_matrix[2127:2442,9]<-'Train2'
gsva_MYC_matrix[2443:2550,9]<-'Test2'
gsva_MYC_matrix[2551:2920,9]<-'Train3'
gsva_MYC_matrix[2921:3084,9]<-'Test3'
#### MYCN-score ####
MYCN_signature_2020<-list(c('ARMC6','DCTPP1','EIF4G1','ELOVL6','FBL','PRMT1','E2F1','RPS19'))
MYCN_signature<-read.csv('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/MYCN-signature.csv',header = T,sep = ',')
MYCN_signature_po<-MYCN_signature[MYCN_signature$Correlation>0,]
MYCN_signature_ne<-MYCN_signature[MYCN_signature$Correlation<0,]
# MYCN_signature<-list(MYCN_signature$Gene)
gsva_MYC_matrix<-as.data.frame(matrix(data=NA,
                                      nrow=3084,ncol=3))
colnames(gsva_MYC_matrix)<-c('Sample','MYCN_score_2020','MYCN_score_PNAS')

gsva_MYC_matrix[1:714,2]<-t(gsva(as.matrix(merge_data_batch_MYCN_amp),MYCN_signature , mx.diff=1,kcdf="Gaussian",
                                  method= "ssgsea",ssgsea.norm=F))

gsva_MYC_matrix[715:1748,2]<-t(gsva(as.matrix(train_1),MYCN_signature_2020 , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[1749:2126,2]<-t(gsva(as.matrix(test_1),MYCN_signature_2020 , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2127:2442,2]<-t(gsva(as.matrix(train_2),MYCN_signature_2020 , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2443:2550,2]<-t(gsva(as.matrix(test_2),MYCN_signature_2020 , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2551:2920,2]<-t(gsva(as.matrix(train_3),MYCN_signature_2020 , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2921:3084,2]<-t(gsva(as.matrix(test_3),MYCN_signature_2020 , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))

gsva_MYC_matrix[1:714,3]<-t(gsva(as.matrix(merge_data_batch_MYCN_amp),list(MYCN_signature_po$Gene) , mx.diff=1,kcdf="Gaussian",
                                  method= "ssgsea",ssgsea.norm=F))

gsva_MYC_matrix[715:1748,3]<-t(gsva(as.matrix(train_1),list(MYCN_signature_po$Gene) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[1749:2126,3]<-t(gsva(as.matrix(test_1),list(MYCN_signature_po$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2127:2442,3]<-t(gsva(as.matrix(train_2),list(MYCN_signature_po$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2443:2550,3]<-t(gsva(as.matrix(test_2),list(MYCN_signature_po$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2551:2920,3]<-t(gsva(as.matrix(train_3),list(MYCN_signature_po$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2921:3084,3]<-t(gsva(as.matrix(test_3),list(MYCN_signature_po$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[1:714,1]<-colnames(merge_data_batch_MYCN_amp)
gsva_MYC_matrix[715:1748,1]<-colnames(train_1)
gsva_MYC_matrix[1749:2126,1]<-colnames(test_1)
gsva_MYC_matrix[2127:2442,1]<-colnames(train_2)
gsva_MYC_matrix[2443:2550,1]<-colnames(test_2)
gsva_MYC_matrix[2551:2920,1]<-colnames(train_3)
gsva_MYC_matrix[2921:3084,1]<-colnames(test_3)

gsva_MYC_matrix[1:714,10]<-t(gsva(as.matrix(merge_data_batch_MYCN_amp),MYCN_signature , mx.diff=1,kcdf="Gaussian",
                                  method= "ssgsea",ssgsea.norm=F))

gsva_MYC_matrix[715:1748,10]<-t(gsva(as.matrix(train_1),list(MYCN_signature$Gene) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[1749:2126,10]<-t(gsva(as.matrix(test_1),list(MYCN_signature$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2127:2442,10]<-t(gsva(as.matrix(train_2),list(MYCN_signature$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2443:2550,10]<-t(gsva(as.matrix(test_2),list(MYCN_signature$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2551:2920,10]<-t(gsva(as.matrix(train_3),list(MYCN_signature$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2921:3084,10]<-t(gsva(as.matrix(test_3),list(MYCN_signature$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))

gsva_MYC_matrix[1:714,12]<-t(gsva(as.matrix(merge_data_batch_MYCN_amp),list(MYCN_signature_ne$Gene) , mx.diff=1,kcdf="Gaussian",
                                 method= "ssgsea",ssgsea.norm=F))

gsva_MYC_matrix[715:1748,12]<-t(gsva(as.matrix(train_1),list(MYCN_signature_ne$Gene) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[1749:2126,12]<-t(gsva(as.matrix(test_1),list(MYCN_signature_ne$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2127:2442,12]<-t(gsva(as.matrix(train_2),list(MYCN_signature_ne$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2443:2550,12]<-t(gsva(as.matrix(test_2),list(MYCN_signature_ne$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2551:2920,12]<-t(gsva(as.matrix(train_3),list(MYCN_signature_ne$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix[2921:3084,12]<-t(gsva(as.matrix(test_3),list(MYCN_signature_ne$Gene) , mx.diff=1,kcdf="Gaussian",
                                      method= "ssgsea",ssgsea.norm=F))

write.csv(gsva_MYC_matrix,'ssGSEA_MYC_MYCN.csv')
gsva_MYC_matrix<-read.csv('GEO+Array/Consensus/Kmean/ssGSEA_MYC_MYCN.csv')
cor.test(gsva_MYC_matrix[,2],log2(gsva_MYC_matrix[,11]),method = 'spearman')
cor.test(gsva_MYC_matrix[,3],log2(gsva_MYC_matrix[,11]),method = 'spearman')
cor.test(gsva_MYC_matrix[,4],log2(gsva_MYC_matrix[,11]),method = 'spearman')

pdf('AURKA-MYC1-spearman.pdf',width = 8,height = 8)# myc1 4500|myc2 3200|myc1+2 4200
par(mar=c(2,4,2,2))
ggscatter(gsva_MYC_matrix, x = "AURKA", y = "MYC1-score",
          size = 1,
          color='Group-subgroup',
          shape='Group-subgroup',
          xlab = expression(paste('AURKA expression')),
          ylab = expression(paste('MYC1',' Score')),
          #add = "reg.line",  # Add regressin line
          #cor.coef = TRUE,
          add.params = list(color = "Group-subgroup"))+ # Customize reg. li
         # conf.int = TRUE, # Add confidence interval
          #cor.coeff.args = list(method = "spearman", label.x =10 , label.sep = "\n",
           #                 label.y = 4500),cor.coef.size = 10)+
  annotate("text", x =7, y=4500, label = substitute(paste(italic("R"), " = 0.70")),size=5,
           color="#F8766D")+
  annotate("text", x =8.2, y=4510, label = substitute(paste(italic("P"), " < 2.2 x ",10^-16)),size=5,
           color="#F8766D")+
  annotate("text", x =7, y=4400, label = substitute(paste(italic("R"), " = 0.71")),size=5,
             color="#7CAE00")+
  annotate("text", x =8.2, y=4410, label = substitute(paste(italic("P"), " < 2.2 x ",10^-16)),size=5,
             color="#7CAE00")+
  annotate("text", x =7, y=4300, label = substitute(paste(italic("R"), " = 0.44")),size=5,
             color="#00BFC4")+
  annotate("text", x =8.2, y=4310, label = substitute(paste(italic("P"), " < 2.2 x ",10^-16)),size=5,
             color="#00BFC4")+
  annotate("text", x =7, y=4200, label = substitute(paste(italic("R"), " = 0.73")),size=5,
             color="#C77CFF")+
  annotate("text", x =8.2, y=4210, label = substitute(paste(italic("P"), " < 2.2 x ",10^-16)),size=5,
             color="#C77CFF")+
  theme(axis.title=element_text(size=20,face="bold"),axis.text=element_text(size=10,face="bold"),
        axis.title.x=element_text(face = 'bold'),
        axis.title.y=element_text(face = 'bold'))
dev.off()
cor.test(gsva_MYC_matrix[,1],log2(gsva_MYC_matrix[,4]),method = 'spearman')

cor.test(gsva_MYC_matrix[1:714,1],log2(gsva_MYC_matrix[1:714,2]),method = 'spearman')
cor.test(gsva_MYC_matrix[715:2126,1],gsva_MYC_matrix[715:2126,11],method = 'spearman')
cor.test(gsva_MYC_matrix[2127:2550,1],gsva_MYC_matrix[2127:2550,11],method = 'spearman')
cor.test(gsva_MYC_matrix[2551:3084,1],gsva_MYC_matrix[2551:3084,2],method = 'spearman')



pdf('AURKA-MYC1+2-spearman-normal+AMP.pdf',width = 8,height = 8)# myc1 4500|myc2 3200|myc1+2 4200
par(mar=c(2,4,2,2))
ggscatter(gsva_MYC_matrix, x = "AURKA", y = "MYC-score",
          size = 1,
          color='black',
          shape='Group',
          xlab = expression(paste('AURKA expression')),
          ylab = expression(paste('MYC1+2',' Score')),
          #add = "reg.line",  # Add regressin line
          #cor.coef = TRUE,
          add.params = list(color = "black"))+ # Customize reg. li
  # conf.int = TRUE, # Add confidence interval
  #cor.coeff.args = list(method = "spearman", label.x =10 , label.sep = "\n",
  #                 label.y = 4500),cor.coef.size = 10)+
  annotate("text", x =7, y=4200, label = substitute(paste(italic("R"), " = 0.83")),size=5,
           color="black")+
  annotate("text", x =8.2, y=4210, label = substitute(paste(italic("P"), " < 2.2 x ",10^-16)),size=5,
           color="black")

dev.off()

my_comparisons_train <- list( c("Subgroup1", "MYCN-AMP"), 
                              c("Subgroup2", "MYCN-AMP"),
                              c("Subgroup3", "MYCN-AMP"),
                              c('Subgroup1','Subgroup2'),
                              c('Subgroup1','Subgroup3'),
                              c('Subgroup2','Subgroup3'))

colnames(Purity_Test)
Purity_Test$Group <- factor(Purity_Test$Group,levels = c("MYCN_AMP", "Subgroup_1", "Subgroup_2",
                                                         "Subgroup_3"))

Purity_Train$Group <- factor(Purity_Train$Group,levels = c("MYCN_AMP", "Subgroup_1", "Subgroup_2",
                                                           "Subgroup_3"))
Purity_All$Group <- factor(Purity_All$Group,levels = c("MYCN_AMP", "Subgroup_1", "Subgroup_2",
                                                       "Subgroup_3"))
229+387

colnames(gsva_MYC_matrix_train)
gsva_MYC_matrix_train<-gsva_MYC_matrix[gsva_MYC_matrix$Train_Test!='Test',]
gsva_MYC_matrix_test<-gsva_MYC_matrix[gsva_MYC_matrix$Train_Test!='Train',]
ggviolin(gsva_MYC_matrix_test, x = "Group.subgroup", y = 'MYCN_po',
            fill = "Group.subgroup",
            add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons_train,size=10)+
  ylab("MYCN score\n") +
  xlab('')+
  #scale_y_continuous(breaks=seq(0, 1, 0.1))+
  #coord_cartesian(ylim = c(-2500,3000))+
  #geom_hline(yintercept = median(gsva_MYC_matrix[Purity_Test$Group=='MYCN_AMP',]$ImmuneScore), 
   #          size=2,linetype = 2)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text( face="bold", size=30,color = 'black',hjust = 0.5,vjust=0),
        axis.text.y = element_text( face="bold", size=30,color = 'black'),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(color="black", size=40, face="bold"),
        plot.margin = unit(c(1,1,1,1),"cm"))
ggsave('Subgroup_AMP_MYCN_test.pdf',width = 16,height = 12) 
rm(gsva_MYC_matrix_train,gsva_MYC_matrix_test)

#### Subgroup2 MYCN and correlated genes ####
### 所有
merge_data_batch_MYCN_amp_train_train1<-data.frame(merge_data_batch_MYCN_amp,
                                                   train_1,train_2,train_3)
dim(merge_data_batch_MYCN_amp_train_train1)

pv<-array(1,c(length(merge_data_batch_MYCN_amp_train_train1[,1]),1)) #  make a empty list for p-value, you should change 13767 to your data value
cor<-array(1,c(length(merge_data_batch_MYCN_amp_train_train1[,1]),1)) # make a empty list for correlation value
for (i in 1:length(merge_data_batch_MYCN_amp_train_train1[,1])) #change 13767
{
  VL<-cor.test(as.numeric(merge_data_batch_MYCN_amp_train_train1[i,]),
               as.numeric(gsva_MYC_matrix_train[,11]),method = "pearson") #### change GRK6 to any gene you want
  pv[i,1]<-VL$p.value    #### extract p-value
  cor[i,1]<-VL$estimate #### extract correlation value
}#log2_cc_ee[i,1]<-log2(mean(c(as.numeric(ee[i,])))/mean(c(as.numeric(cc[i,]))))}
#### 单独subgroup
gsva_MYC_matrix_train_sub<-gsva_MYC_matrix_train[gsva_MYC_matrix_train$Group.subgroup=='MYCN-AMP',]
pv<-array(1,c(length(merge_data_batch_MYCN_amp[,1]),1)) #  make a empty list for p-value, you should change 13767 to your data value
cor<-array(1,c(length(merge_data_batch_MYCN_amp[,1]),1)) # make a empty list for correlation value
for (i in 1:length(merge_data_batch_MYCN_amp[,1])) #change 13767
{
  VL<-cor.test(as.numeric(merge_data_batch_MYCN_amp[i,]),
               as.numeric(gsva_MYC_matrix_train_sub[,11]),method = "pearson") #### change GRK6 to any gene you want
  pv[i,1]<-VL$p.value    #### extract p-value
  cor[i,1]<-VL$estimate #### extract correlation value
}#log2_cc_ee[i,1]<-log2(mean(c(as.numeric(ee[i,])))/mean(c(as.numeric(cc[i,]))))}

fdr<-p.adjust(pv, method = 'BH', n = length(pv)) ## adjust p-value
fdr<-as.data.frame(fdr)

gsva_MYC_matrix_train_sub<-data.frame(rownames(merge_data_batch_MYCN_amp_train_train1),
                                  pv,fdr,cor)#,log2_cc_ee)   # create a new data frame with genenames, p-value and correlation value

write.table(gsva_MYC_matrix_train_sub,'AMP_Train_MYCN_correlation.txt',sep = '\t')
rm(gsva_MYC_matrix_train_sub,pv,VL,cor,fdr)
### Density plot of correlation coefficient
AMP_Train_MYCN_correlation<-read.table('GEO+Array/Consensus/Kmean/Subgroup2_MYCN/AMP_Train_MYCN_correlation.txt',
                                       sep = '\t',header = T)

dim(AMP_Train_MYCN_correlation)
AMP_Train_MYCN_correlation_density<-as.data.frame(matrix(data=NA,
                                                         nrow=11583*4,ncol=2))#四个组coefficient的matrix
11583*4
AMP_Train_MYCN_correlation_density[1:11583,1]<-AMP_Train_MYCN_correlation$COR_Train1
AMP_Train_MYCN_correlation_density[1:11583,2]<-'Subgroup1'
AMP_Train_MYCN_correlation_density[11584:23166,1]<-AMP_Train_MYCN_correlation$COR_Train2
AMP_Train_MYCN_correlation_density[11584:23166,2]<-'Subgroup2'
AMP_Train_MYCN_correlation_density[23167:34749,1]<-AMP_Train_MYCN_correlation$COR_Train3
AMP_Train_MYCN_correlation_density[23167:34749,2]<-'Subgroup3'
AMP_Train_MYCN_correlation_density[34750:46332,1]<-AMP_Train_MYCN_correlation$COR_AMP
AMP_Train_MYCN_correlation_density[34750:46332,2]<-'MYCN-AMP'
colnames(AMP_Train_MYCN_correlation_density)<-c('Correlation','Group')
p1 <- ggplot(data=AMP_Train_MYCN_correlation_density, aes(x=Correlation, group=Group, fill=Group))+
  xlab("Pearson's correlation") + ylab("Denisty") +
  geom_density(adjust=2,alpha=.6) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  theme_classic()
p1<-p1+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.margin =unit(c(1,1,1,1), "cm"))
p1
ggsave('AMP_Train_MYCN_correlation.pdf',width = 15,height = 6) #保存图片
rm(AMP_Train_MYCN_correlation_density,p1)
AMP_Train_MYCN_correlation_Train1<-subset(AMP_Train_MYCN_correlation,
                                          FDR_Train1 < 0.05 & COR_Train1 > 0.5)
AMP_Train_MYCN_correlation_Train2<-subset(AMP_Train_MYCN_correlation,
                                          FDR_Train2 < 0.05 & COR_Train2 > 0.5)
AMP_Train_MYCN_correlation_Train3<-subset(AMP_Train_MYCN_correlation,
                                          FDR_Train3 < 0.05 & COR_Train3 > 0.5)
AMP_Train_MYCN_correlation_AMP<-subset(AMP_Train_MYCN_correlation,
                                       FDR_AMP < 0.05 & COR_AMP > 0.5)
BiocManager::install('VennDiagram')
library(VennDiagram)

venn.diagram(list(MYCN_AMP=as.factor(AMP_Train_MYCN_correlation_AMP$Gene),
                  Subgroup1=as.factor(AMP_Train_MYCN_correlation_Train1$Gene),
                  Subgroup2=as.factor(AMP_Train_MYCN_correlation_Train2$Gene),
                  Subgroup3=as.factor(AMP_Train_MYCN_correlation_Train3$Gene)),
             resolution = 300, imagetype = "tiff", 
             col='transparent',
             fill=c("lightblue3", "mediumpurple2", "red3",'goldenrod1'),
             main="Differentially expressed proteins",
             main.cex = 3, main.fontface = 2,alpha = 0.50, cex=2.5, cat.cex=2.5,
             margin=0.2, #改变图片大小
             cat.dist = c(0.24, 0.24, 0.1,0.1),##注释距离圆的距离
             # cex for each area label
             #cat.cex for category name
             #cat.col = c("darkblue", "darkgreen", "orange"),
             #print.mode = 'percent',
             #category.names = c('RAS(75)',expression(paste('TGF-',beta,"(35)")),expression(paste('RAS+TGF-',beta,"(47)"))),
             filename = "AMP_Train_MYCN_score_positive.tiff")

intersect(AMP_Train_MYCN_correlation_AMP$Gene,MYCN_signature_po$Gene)
####Subgroup2与1/3的DEG和MYCN-pos-correlated gene的overlap
IPA_subgroup2DEG<-read.table('GEO+Array/Consensus/Kmean/IPA_Subgroup1.txt',sep = '\t',header = T)
IPA_subgroup2DEG<-IPA_subgroup2DEG[,c(1,7,11)]
IPA_subgroup2DEG_pos<-subset(IPA_subgroup2DEG,adj.P.Valsub2< 0.05 & logFCsub2>1)
IPA_subgroup2DEG_neg<-subset(IPA_subgroup2DEG,adj.P.Valsub2< 0.05 & logFCsub2 < -1)
AMP_Train_MYCN_correlation_Train2_pos<-subset(AMP_Train_MYCN_correlation,
                                          FDR_Train2 < 0.05 & COR_Train2 > 0.3)
AMP_Train_MYCN_correlation_Train2_neg<-subset(AMP_Train_MYCN_correlation,
                                          FDR_Train2 < 0.05 & COR_Train2 < -0.3)

venn.diagram(list(Increased=as.factor(IPA_subgroup2DEG_pos$Gene),
                  Positive=as.factor(AMP_Train_MYCN_correlation_Train2_pos$Gene),
                  Decreased=as.factor(IPA_subgroup2DEG_neg$Gene),
                  Negative=as.factor(AMP_Train_MYCN_correlation_Train2_neg$Gene)),
             resolution = 300, imagetype = "tiff", 
             col='transparent',
             fill=c("lightblue3", "mediumpurple2", "red3",'goldenrod1'),
             main="Differentially expressed proteins",
             main.cex = 3, main.fontface = 2,alpha = 0.50, cex=2.5, cat.cex=2.5,
             margin=0.2, #改变图片大小
             cat.dist = c(0.24, 0.24, 0.1,0.1),##注释距离圆的距离
             # cex for each area label
             #cat.cex for category name
             #cat.col = c("darkblue", "darkgreen", "orange"),
             #print.mode = 'percent',
             #category.names = c('RAS(75)',expression(paste('TGF-',beta,"(35)")),expression(paste('RAS+TGF-',beta,"(47)"))),
             filename = "Subgroup2_DEG_Correlation.tiff")
Subgroup2_Up_Pos<-intersect(IPA_subgroup2DEG_pos$Gene,AMP_Train_MYCN_correlation_Train2_pos$Gene)

Subgroup2_Up_Pos_removeMYCNSignature<-setdiff(Subgroup2_Up_Pos,MYCN_signature$Gene)
Subgroup2_Up_Pos-MYCN_signature_po$Gene
write.csv(Subgroup2_Up_Pos,'Subgroup2_Up_Pos_all.csv')
write.csv(Subgroup2_Up_Pos_removeMYCNSignature,'Subgroup2_Up_Pos_removeMYCNSignature.csv')
### Subgroup2 negative gene 
Subgroup2_Up_Neg<-intersect(IPA_subgroup2DEG_neg$Gene,AMP_Train_MYCN_correlation_Train2_neg$Gene)

Subgroup2_Up_Neg_removeMYCNSignature<-setdiff(Subgroup2_Up_Neg,MYCN_signature$Gene)  ###无显著PPI网络

AMP_Train_MYCN_correlation_Train2_neg

write.table(AMP_Train_MYCN_correlation_Train2_neg$Gene,'AMP_Train_MYCN_correlation_Train2_neg.txt')
AMP_train_neg<-read.xlsx('/Users/yihuawang/Downloads/all.tfeawagc1/metascape_result.xlsx')
sort(table(AMP_train_neg$Gene_Location))

mean(train_2_sil['AURKA',])
mean(train_1['EZH2',])
mean(train_3['EZH2',])
mean(as.matrix(merge_data_batch_MYCN_amp)['EZH2',])
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

gsva_Immunepathway_matrix[1:714,51]<-colnames(merge_data_batch_MYCN_amp)
gsva_Immunepathway_matrix[715:1748,51]<-colnames(train_1)
gsva_Immunepathway_matrix[1749:2126,51]<-colnames(test_1)
gsva_Immunepathway_matrix[2127:2442,51]<-colnames(train_2)
gsva_Immunepathway_matrix[2443:2550,51]<-colnames(test_2)
gsva_Immunepathway_matrix[2551:2920,51]<-colnames(train_3)
gsva_Immunepathway_matrix[2921:3084,51]<-colnames(test_3)
colnames(gsva_Immunepathway_matrix)[50:51]<-c('Train_Test','Sample')
write.csv(gsva_Immunepathway_matrix,'gsva_Immunepathway_matrix.csv')
gsva_Immunepathway_matrix_train<-gsva_Immunepathway_matrix[gsva_Immunepathway_matrix$Train_Test!='Test',]
for (i in colnames(gsva_Immunepathway_matrix)[1:46]){
  print(i)
  ggviolin(gsva_Immunepathway_matrix_train, x = "Group-subgroup", y = i,
           fill = "Group-subgroup",
           add = "boxplot", add.params = list(fill = "white"))+
    stat_compare_means(comparisons = my_comparisons_train,size=10)+
    ylab(paste0(i,'\n')) +
    xlab('')+
    #scale_y_continuous(breaks=seq(0, 1, 0.1))+
    #coord_cartesian(ylim = c(-2500,3000))+
    #geom_hline(yintercept = median(gsva_MYC_matrix[Purity_Test$Group=='MYCN_AMP',]$ImmuneScore), 
    #          size=2,linetype = 2)+
    theme_classic()+
    theme(legend.position = "none",
          axis.text.x = element_text( face="bold", size=30,color = 'black',hjust = 0.5,vjust=0),
          axis.text.y = element_text( face="bold", size=30,color = 'black'),
          axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
          axis.title.y = element_text(color="black", size=40, face="bold"),
          plot.margin = unit(c(1,1,1,1),"cm"))
  ggsave(paste0(i,'_ImmunePathway_AMP_MYCN_train.pdf'),width = 16,height = 12)
}

rownames(gsva_Immunepathway_matrix_train)<-gsva_Immunepathway_matrix_train$Sample
gsva_Immunepathway_matrix_train<-gsva_Immunepathway_matrix_train[,1:46]
gsva_Immunepathway_matrix_train<-t(gsva_Immunepathway_matrix_train)
n=t(scale(t(gsva_Immunepathway_matrix_train)))
n[n>3]=3
n[n< -3]= -3
n[1:4,1:4]
ac<-data.frame(c(rep('MYCN-AMP',714),rep('Subgroup1',1034),rep('Subgroup2',316), # AMP+Train
                rep('Subgroup3',370)))

ac<-data.frame(c(rep('Subgroup1',1034),rep('Subgroup2',316), # AMP+Train
                 rep('Subgroup3',370)))

ac<-data.frame(c(rep('Subgroup3',370)))

714+1034+316
colnames(ac)<-'Molecular subtype'
#ac$`Molecular subtype` <- factor(ac$`Molecular subtype`, levels = c("MYCN-AMP", "Subgroup1",'Subgroup2','Subgroup3'))
rownames(ac)=colnames(n)[715:2434]
rownames(ac)=colnames(n)[2065:2434]

ac$Blank<-'1'
ac<-as.matrix(t(ac))
dim(n_order)
dim(ac)
set.seed(123)
pdf("heatmap_Immunepathway_order.pdf", width = 10, height = 10)
ha = HeatmapAnnotation(`Molecular subtype` =ac[1,],
                       col = list(`Molecular subtype` = c('MYCN-AMP' = "#F8766D", 'Subgroup1' = "#7CAE00",
                                                          'Subgroup2'='#00BFC4','Subgroup3'='#C77CFF')),
                       annotation_legend_param = list(
                         `Molecular subtype` = list(title = "Molecular subtype")))
draw(Heatmap(as.matrix(n_order), name = "Immune pathway", 
             column_km = 1,column_km_repeats = 100,
             clustering_distance_columns = "euclidean",
             column_order = F,
             show_row_names = T,show_column_names  = FALSE,
             top_annotation = ha,col=colorbar,show_heatmap_legend = T,
             show_parent_dend_line = FALSE),padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()
heatmap_subgroup1_immune<-Heatmap(n[,2065:2434], name = "Immune pathway", 
                                  column_km = 1,column_km_repeats = 1000,
                                  clustering_distance_columns = "euclidean",
                                  show_row_names = FALSE,show_column_names  = FALSE,
                                  top_annotation = ha,col=colorbar,show_heatmap_legend = T,
                                  show_parent_dend_line = FALSE)

heatmap_subgroup1_immune = draw(heatmap_subgroup1_immune)
row_order(heatmap_subgroup1_immune)

n_AMP<-n[,2065:2434]

train_AMP_immune<-n_AMP[,column_order(heatmap_subgroup1_immune)]
train_1_immune<-n_AMP[,column_order(heatmap_subgroup1_immune)]
train_2_immune<-n_AMP[,column_order(heatmap_subgroup1_immune)]
train_3_immune<-n_AMP[,column_order(heatmap_subgroup1_immune)]
n_order<-data.frame(train_AMP_immune,train_1_immune,
                    train_2_immune,train_3_immune)

rownames(ac)=colnames(n_order)
ac$Blank<-'1'
ac<-as.matrix(t(ac))

pdf("heatmap_Immunepathway_order.pdf", width = 15, height = 10)
ha = HeatmapAnnotation(`Molecular subtype` =ac[1,],
                       col = list(`Molecular subtype` = c('MYCN-AMP' = "#F8766D", 'Subgroup1' = "#7CAE00",
                                                          'Subgroup2'='#00BFC4','Subgroup3'='#C77CFF')),
                       annotation_legend_param = list(
                         `Molecular subtype` = list(title = "Molecular subtype")))
draw(Heatmap(as.matrix(n_order), name = "Immune pathway", 
             column_km = 1,column_km_repeats = 100,
             clustering_distance_columns = "euclidean",
             cluster_columns = FALSE,
             show_row_names = T,show_column_names  = FALSE,
             top_annotation = ha,col=colorbar,show_heatmap_legend = T,
             show_parent_dend_line = FALSE),padding = unit(c(2, 2, 2, 2), "cm"))
dev.off()
train_1['PRF1',1]
rm(train_AMP_immune,train_1_immune,train_2_immune,train_3_immune)
#### correlation between Immune and hypoxia TGFB and KRAS 	PI3K-Akt MAPK，RAS signature ####
Immune_genelist<-list(c('LCP2',	'LSP1',	'FYB',	'PLEK'	,'HCK',	'IL10RA',	'LILRB1',	'NCKAP1L'	,'LAIR1',	'NCF2',
                        'CYBB','PTPRC',	'IL7R',	'LAPTM5',	'CD53',	'EVI2B',	'SLA',	'ITGB2',	'GIMAP4',	'MYO1F',	
                        'HCLS1',	'MNDA',	'IL2RG',	'CD48',	'AOAH',	'CCL5',	'LTB',	'GMFG',	'GIMAP6',	'GZMK',	'LST1',
                        'GPR65',	'LILRB2',	'WIPF1',	'CD37',	'BIN2',	'FCER1G',	'IKZF1',	'TYROBP',	'FGL2',	'FLI1',
                        'IRF8',	'ARHGAP15',	'SH2B3',	'TNFRSF1B',	'DOCK2',	'CD2',	'ARHGEF6',	'CORO1A',	'LY96',
                        'LYZ',	'ITGAL',	'TNFAIP3',	'RNASE6',	'TGFB1',	'PSTPIP1',	'CST7',	'RGS1',	'FGR',	'SELL',
                        'MICAL1',	'TRAF3IP3',	'ITGA4',	'MAFB',	'ARHGDIB',	'IL4R',	'RHOH',	'HLA-DPA1',	'NKG7',	'NCF4',
                        'LPXN',	'ITK',	'SELPLG',	'HLA-DPB1',	'CD3D',	'CD300A',	'IL2RB',	'ADCY7',	'PTGER4',	'SRGN',	
                        'CD247',	'CCR7',	'MSN',	'ALOX5AP',	'PTGER2',	'RAC2',	'GBP2',	'VAV1',	'CLEC2B',	'P2RY14',	
                        'NFKBIA',	'S100A9',	'IFI30',	'MFSD1',	'RASSF2',	'TPP1',	'RHOG',	'CLEC4A',	'GZMB',	'PVRIG',
                        'S100A8',	'CASP1',	'BCL2A1',	'HLA-E',	'KLRB1',	'GNLY',	'RAB27A',	'IL18RAP',	'TPST2',	
                        'EMP3',	'GMIP',	'LCK',	'IL32',	'PTPRCAP',	'LGALS9',	'CCDC69',	'SAMHD1',	'TAP1',	'GBP1',
                        'CTSS',	'GZMH',	'ADAM8',	'GLRX',	'PRF1',	'CD69',	'HLA-B',	'HLA-DMA',	'CD74',
                        'KLRK1',	'PTPRE',	'HLA-DRA',	'VNN2',	'TCIRG1',	'RABGAP1L',	'CSTA',	'ZAP70',
                        'HLA-F',	'HLA-G',	'CD52',	'CD302',	'CD27'))

BiocManager::install('msigdbr')
library(msigdbr)
?msigdbr
msigdbr_HS<-msigdbr(species = "Homo sapiens")
table(msigdbr_HS$gs_cat)
length(unique(msigdbr_HS_h$gs_name))

## H, C2, and C5 
unique(msigdbr_HS$gs_cat)
msigdbr_HS_h<-subset(msigdbr_HS,msigdbr_HS$gs_cat=='H')

gsva_MYC_matrix<-as.data.frame(matrix(nrow = 370,ncol = 2,data = NA))
gsva_MYC_matrix$V1<-colnames(train_3)
gsva_MYC_matrix$V2<-t(gsva(as.matrix(train_3),Immune_genelist , mx.diff=1,kcdf="Gaussian",
                           method= "ssgsea",ssgsea.norm=F))
colnames(gsva_MYC_matrix)<-c('Sample','Immune')
gsva_Subgroup3_Immune_Oncogenic_matrix<-as.data.frame(matrix(nrow = 100,ncol = 4,data = NA))

list_a<-1
for (j in c('CP:KEGG','CP:REACTOME','CP:BIOCARTA','BP')){
  msigdbr_HS_h<-subset(msigdbr_HS,msigdbr_HS$gs_subcat==j)
  for (i in unique(msigdbr_HS_h$gs_name)){
    list_msigdbr_HS_h<-list(msigdbr_HS_h[msigdbr_HS_h$gs_name==i,]$gene_symbol)
    gsva_MYC_matrix_i<-t(gsva(as.matrix(train_3),list_msigdbr_HS_h , mx.diff=1,kcdf="Gaussian",
                              method= "ssgsea",ssgsea.norm=F))
    
    cortest<-cor.test(gsva_MYC_matrix_i,gsva_MYC_matrix$Immune,method = 'spearman')
    if (cortest$p.value < 0.05){
      gsva_MYC_matrix$`i`<-gsva_MYC_matrix_i
      gsva_Subgroup3_Immune_Oncogenic_matrix[list_a,]<-c(i,cortest$estimate,cortest$p.value,j)
      list_a<-list_a+1
    }
  }
}
unique(msigdbr_HS$gs_cat)
for (j in c('C6','H')){
  msigdbr_HS_h<-subset(msigdbr_HS,msigdbr_HS$gs_cat==j)
  for (i in unique(msigdbr_HS_h$gs_name)){
    list_msigdbr_HS_h<-list(msigdbr_HS_h[msigdbr_HS_h$gs_name==i,]$gene_symbol)
    gsva_MYC_matrix_i<-t(gsva(as.matrix(train_3),list_msigdbr_HS_h , mx.diff=1,kcdf="Gaussian",
                              method= "ssgsea",ssgsea.norm=F))
    
    cortest<-cor.test(gsva_MYC_matrix_i,gsva_MYC_matrix$Immune,method = 'spearman')
    if (cortest$p.value < 0.05){
      gsva_MYC_matrix$`i`<-gsva_MYC_matrix_i
      gsva_Subgroup3_Immune_Oncogenic_matrix[list_a,]<-c(i,cortest$estimate,cortest$p.value,j)
      list_a<-list_a+1
    }
  }
}
for (j in c('BP')){
  msigdbr_HS_h<-subset(msigdbr_HS,msigdbr_HS$gs_subcat==j)
  for (i in unique(msigdbr_HS_h$gs_name)){
    list_msigdbr_HS_h<-list(msigdbr_HS_h[msigdbr_HS_h$gs_name==i,]$gene_symbol)
    gsva_MYC_matrix_i<-t(gsva(as.matrix(train_3),list_msigdbr_HS_h , mx.diff=1,kcdf="Gaussian",
                              method= "ssgsea",ssgsea.norm=F))
    
    cortest<-cor.test(gsva_MYC_matrix_i,gsva_MYC_matrix$Immune,method = 'spearman')
    print(cortest$p.value)
    if (cortest$p.value < 0.05){
      gsva_MYC_matrix$`i`<-gsva_MYC_matrix_i
      gsva_Subgroup3_Immune_Oncogenic_matrix[list_a,]<-c(i,cortest$estimate,cortest$p.value,j)
      list_a<-list_a+1
    }
  }
}


write.csv(gsva_Subgroup3_Immune_Oncogenic_matrix,'
          gsva_Subgroup3_Immune_Oncogenic_matrix.csv')

gsva_MYC_matrix_i<-t(gsva(as.matrix(train_3),Hypoxia_list , mx.diff=1,kcdf="Gaussian",
                          method= "ssgsea",ssgsea.norm=F))

cortest<-cor.test(gsva_MYC_matrix_i,gsva_MYC_matrix$Immune,method = 'spearman')
rm(gsva_MYC_matrix1)

gsva_MYC_matrix<-na.omit(gsva_MYC_matrix)
Purity_All<-read.csv('/Users/lefan/Rcode/Neuroblastoma/Hu/GEO+Array/Tumorpurity.csv',header = T,sep = ',')
rownames(Purity_All)<-Purity_All$Sample
Purity_All<-Purity_All[,-5] # for score
Purity_All<-na.omit(Purity_All) # for purity


Purity_subgroup3<-Purity_All[colnames(train_3),]
Purity_subgroup3<-na.omit(Purity_subgroup3)

Hypoxia_list<-list(c('NDRG1','ENO1','VEGFA','MRPS17','TPI1','CDKN3',
                     'MIF','LDHA','ALDOA','TUBB6','PGAM1','SLC2A1','P4HA1','ACOT7','ADM'))
gsva_MYC_matrix<-t(gsva(as.matrix(train_3),Hypoxia_list , mx.diff=1,kcdf="Gaussian",
                                 method= "ssgsea",ssgsea.norm=F))
gsva_MYC_matrix<-as.data.frame(gsva_MYC_matrix)
gsva_MYC_matrix$Sample<-rownames(gsva_MYC_matrix)
colnames(gsva_MYC_matrix)[1]<-'Hypoxia'
gsva_MYC_matrix<-merge(gsva_MYC_matrix,Purity_subgroup3)
cor.test(gsva_MYC_matrix$Hypoxia,gsva_MYC_matrix$ImmuneScore,method = 'spearman')
plot(gsva_MYC_matrix$Hypoxia,gsva_MYC_matrix$ImmuneScore)
#### MHC-1 socre based on mean expression of HLA-A, HLA-B, HLA-C, TAP1, TAP2, NLRC5, PSMB9, PSMB8, and B2M
merge_data_batch_MYCN_amp_train_train1<-data.frame(merge_data_batch_MYCN_amp,train_1_sil,train_2_sil,train_3_sil)
dim(merge_data_batch_MYCN_amp_train_train1)
rownames(gsva_Immunepathway_matrix_train)
MHC_CYT<-data.frame(matrix(data=NA,ncol=4,nrow=2434))
colnames(MHC_CYT)<-c('Sample','Group','MHC-1','CYT')
MHC_CYT$Sample<-colnames(merge_data_batch_MYCN_amp_train_train1)
MHC_CYT[1:714,2]<-'MYCN-AMP'
MHC_CYT[715:1748,2]<-'Subgroup1'
MHC_CYT[1749:2064,2]<-'Subgroup2'
MHC_CYT[2065:2434,2]<-'Subgroup3'
train_1['ALDOA',1]
### 无HLA-A,HLA-B,HLA-C
merge_data_batch_MYCN_amp_train_train1_MHC<-merge_data_batch_MYCN_amp_train_train1[c('HLA-A', 'HLA-B', 
                                                                                     'HLA-C', 'TAP1', 'TAP2', 
                                                                                     'NLRC5', 'PSMB9',
                                                                                     'PSMB8','B2M'),]
merge_data_batch_MYCN_amp_train_train1_MHC[,1]
MHC_CYT$`MHC-1`<-apply(merge_data_batch_MYCN_amp_train_train1_MHC[4:9,],2,mean)
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
write.csv(MHC_CYT,'MHC_CYT.csv')
#### Cox regression ####

survival_data_AMP_Cox<-merge(survival_data_AMP,Clinical_test_normal,by='Sample')
survival_data_AMP_Cox<-merge(survival_data_AMP_Cox,gsva_MYC_matrix,by='Sample')
survival_data_AMP_Cox<-merge(survival_data_AMP,Clinical_test_normal,by='Sample')
dim(survival_data_AMP_Cox)

write.csv(survival_data_AMP_Cox,'survival_data_AMP_Cox_AURKA_MYC_all.csv')# <18:1,>18:2| female:1,male:2 | lowrisk:1,high:2|subtype1:1,2:3,3:2| inss 4s:1,1:2,2:3,3:4,4:5
survival_data_AMP_Cox<-read.csv('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/Survival/Cox/survival_data_AMP_Cox_AURKA_MYC_all.csv',row.names = 1,header = T,stringsAsFactors = F)
colnames(survival_data_AMP_Cox)
single_line<-Surv(time =survival_data_AMP_Cox$Day_Survival,event =survival_data_AMP_Cox$Event)
colnames(survival_data_AMP_Cox)[10:12]<-c('MYC_TARGETS_V1','MYC_TARGETS_V2','MYC_TARGETS_V1+V2')
survival_data_AMP_Cox$single_line<- with(survival_data_AMP_Cox,single_line)
Gcox<-coxph(single_line~Subgroup+Age+INSS.Stage+Gender+High.risk+AURKA+
              MYC_TARGETS_V1+MYC_TARGETS_V2+`MYC_TARGETS_V1+V2` ,data = survival_data_AMP_Cox)
ggforest(Gcox,data = survival_data_AMP_Cox,cpositions = c(0.02,0.22,0.4),main = 'Hazard ratio',
         fontsize =1.8,refLabel = 'reference',noDigits = 2)
dev.new()

summary(Gcox)
table(survival_data_AMP_Cox$High.risk)
#### MYCN score-CORRELATED-GENES ####
library(plyr)
library(survival)
colnames(survival_data_AMP_Cox)
sub2_up_pos_removeMYCN<-read.csv('sub2-up-pos_removeMYCN.csv',header = T,sep = ',')# 获取mycn相关的基因 
merge_data_batch_MYCN_amp_train_train1<-data.frame(train_1,train_2,train_3) ### 整合train data
##提取 mycn相关基因在train data中的表达量
sub2_up_pos_removeMYCN_matrix<-merge_data_batch_MYCN_amp_train_train1[sub2_up_pos_removeMYCN$Up.DEG.in.Subgroup2...Correlation.0.3..FDR..0.05..remove.MYCN.siganture,]
sub2_up_pos_removeMYCN_matrix_t<-as.data.frame(t(sub2_up_pos_removeMYCN_matrix))
sub2_up_pos_removeMYCN_matrix_t<-apply(sub2_up_pos_removeMYCN_matrix_t, 2, FUN1)
FUN1
for (j in 1:116){
  col_mean<-mean(sub2_up_pos_removeMYCN_matrix_t[,j])
  for (i in 1: length(sub2_up_pos_removeMYCN_matrix_t[,j])){
    if (sub2_up_pos_removeMYCN_matrix_t[i,j] < col_mean ){
      sub2_up_pos_removeMYCN_matrix_t[i,j]<- 1
    } else {
      sub2_up_pos_removeMYCN_matrix_t[i,j]<-2
    }
  }
}
survival_data<-survival_data[c(rownames(Clinical_test_normal),
                               colnames(merge_data_batch_MYCN_amp_train_train1)),]

survival_data<-na.omit(survival_data)
survival_data$Group<-NA
survival_data[colnames(train_1),4]<-'Subgroup1' ## 6Group所在的列
survival_data[colnames(train_2),4]<-'Subgroup2' ## 6是根据Group所在的列
survival_data[colnames(train_3),4]<-'Subgroup3' ## 6是根据Group所在的列
survival_data_AMP_Cox<-merge(survival_data_AMP,Clinical_test_normal,by='Sample')
sub2_up_pos_removeMYCN_matrix_t$Sample<-rownames(sub2_up_pos_removeMYCN_matrix_t)
survival_data_AMP_Cox<-merge(survival_data_AMP_Cox,sub2_up_pos_removeMYCN_matrix_t,by='Sample')
survival_data_AMP_Cox_surv<-survival_data_AMP_Cox[,c(1,2,3,10:125)]
Basurv<-Surv(time=survival_data_AMP_Cox_surv$Day_Survival,event = survival_data_AMP_Cox_surv$Event)
survival_data_AMP_Cox_surv$Basurv<-with(survival_data_AMP_Cox_surv,Basurv)

UniCox<-function(x){
  FML<-as.formula(paste0('Basurv~',x))
  Gcox<-coxph(FML,data=survival_data_AMP_Cox)
  Gsum<-summary(Gcox)
  
  HR<-round(Gsum$coefficients[,2],3)
  Pvalue<-Gsum$coefficients[,5]
  CI<-paste0(round(Gsum$conf.int[,3:4],2),collapse = '-')
  Unicox<-data.frame('Characteristics'=x,
                     'Hazard Ratio'=HR,
                     'CI95'=CI,
                     'P value'=Pvalue)
  return(Unicox)
}
VarNames<-colnames(survival_data_AMP_Cox_surv)[c(4:119)]
UniVar<-lapply(VarNames, UniCox)
UniVar<-ldply(UniVar,data.frame)
getwd()
write.csv(UniVar,'Subgrooup2_Unicox.csv')

Basurv<-Surv(time=survival_data_AMP_Cox$Day_Survival,event = survival_data_AMP_Cox$Event)
survival_data_AMP_Cox$Basurv<-with(survival_data_AMP_Cox,Basurv)
Gcox<-coxph(Basurv~Subgroup+High.risk+AURKA ,data = survival_data_AMP_Cox)

library( survminer)
pdf('gcox.pdf',height = 10,width = 10)
ggforest(Gcox,data = survival_data_AMP_Cox,cpositions = c(0.02,0.22,0.4),main = 'Hazard ratio',
         fontsize =1.8,refLabel = 'reference',noDigits = 2)
dev.off()
# AURKA 9.291258 |MYC1 3866.276| MYC2 1813.611 | MYC1+2 3415.986 | MYCN-2020 3040.509 | MYCN-PNAS 1402.215

median(as.matrix(gsva_MYC_matrix[715:3084,4]))
median(as.matrix(gsva_MYC_matrix[715:3084,3]))
survival_data_AMP_Cox_AURKA_MYCN_all<-read.csv('Survival/Cox/survival_data_AMP_Cox_AURKA_MYC_all.csv',
                                               header = T,sep = ',')
survival_data_AMP_Cox_AURKA_MYCN_all<-merge(survival_data_AMP_Cox_AURKA_MYCN_all,
                                            gsva_MYC_matrix,by='Sample')
write.csv(survival_data_AMP_Cox_AURKA_MYCN_all,'survival_data_AMP_Cox_AURKA_MYCN_all.csv')
survival_data_AMP_Cox_AURKA_MYCN_all<-read.csv('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/survival_data_AMP_Cox_AURKA_MYCN_all copy.csv',
                                               header = T,sep = ',')
class(survival_data_AMP_Cox_AURKA_MYCN_all[7,9])
gsva_MYC_matrix

survival_data_AMP_Cox_AURKA_MYCN_Risk<-survival_data_AMP_Cox_AURKA_MYCN_all
survival_data_AMP_Cox_AURKA_MYCN_Risk<-merge(survival_data_AMP_Cox_AURKA_MYCN_Risk,gsva_MYC_matrix[,c(1,2,11)],by='Sample')
survival_data_AMP_Cox_AURKA_MYCN_Risk<-survival_data_AMP_Cox_AURKA_MYCN_Risk[survival_data_AMP_Cox_AURKA_MYCN_Risk$High_risk!='Not defined',]
survival_data_AMP_Cox_AURKA_MYCN_Risk<-survival_data_AMP_Cox_AURKA_MYCN_Risk[survival_data_AMP_Cox_AURKA_MYCN_Risk$Subgroup=='Subgroup1',]
### 只包含 low/high risk
survival_data_AMP_Cox_AURKA_MYCN_Risk<-survival_data_AMP_Cox_AURKA_MYCN_all[survival_data_AMP_Cox_AURKA_MYCN_all$High_risk!='Not defined',]
survival_data_AMP_Cox_AURKA_MYCN_Risk<-survival_data_AMP_Cox_AURKA_MYCN_all[survival_data_AMP_Cox_AURKA_MYCN_all$Age!='Not defined',]
survival_data_AMP_Cox_AURKA_MYCN_Risk<-survival_data_AMP_Cox_AURKA_MYCN_all[survival_data_AMP_Cox_AURKA_MYCN_all$Subgroup=='Subgroup2',]
survival_data_AMP_Cox_AURKA_MYCN_Risk$High_risk <- factor(survival_data_AMP_Cox_AURKA_MYCN_Risk$High_risk,levels = c('Low risk','High risk'))
survival_data_AMP_Cox_AURKA_MYCN_Risk$Gender <- factor(survival_data_AMP_Cox_AURKA_MYCN_Risk$Gender,levels = c('Male','Female'))

survival_data_AMP_Cox_AURKA_MYCN_Risk$INSS_Stage <- factor(survival_data_AMP_Cox_AURKA_MYCN_Risk$INSS_Stage,
                                                           levels = c('Stage2','Stage3','Stage4','Stage4S','Stage1'))


Basurv_Risk<-Surv(time=survival_data_AMP_Cox_AURKA_MYCN_Risk$Day_Survival,event = survival_data_AMP_Cox_AURKA_MYCN_Risk$Event)
# 单因分析
table(survival_data_AMP_Cox_AURKA_MYCN_Risk$AURKA.y)
colnames(survival_data_AMP_Cox_AURKA_MYCN_Risk)
Gcox_Risk<-coxph(Basurv_Risk~Age+INSS_Stage+MYCN_score_PNAS,data = survival_data_AMP_Cox_AURKA_MYCN_Risk)
Gcox_Risk<-coxph(Basurv_Risk~AURKA.y,data = survival_data_AMP_Cox_AURKA_MYCN_Risk)
Gcox_Risk<-coxph(Basurv_Risk~Age,data = survival_data_AMP_Cox_AURKA_MYCN_Risk)
summary(Gcox_Risk)
dim(Gcox_Risk)
Gcox_Risk$coefficients
dev.new()
pdf('Gcox_Risk_MYCN_score_PNAS_317sample.pdf',width=15,height = 6)
ggforest(Gcox_Risk,data = survival_data_AMP_Cox_AURKA_MYCN_Risk,cpositions = c(0.02,0.22,0.4),main = 'Hazard ratio',
         fontsize =1.8,refLabel = 'reference',noDigits = 2)
dev.off()
?ggforest
HR <-signif(x$coef[2], digits=2);#exp(beta)
HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)

Basurv<-Surv(time=survival_data_AMP_Cox_AURKA_MYCN_all$Day_Survival,event = survival_data_AMP_Cox_AURKA_MYCN_all$Event)
#survival_data_AMP_Cox_AURKA_MYCN_all$Basurv<-with(survival_data_AMP_Cox_AURKA_MYCN_all,Basurv)
survival_data_AMP_Cox_AURKA_MYCN_all$High_risk <- factor(survival_data_AMP_Cox_AURKA_MYCN_all$High_risk,levels = c('Not defined','Low risk','High risk'))
survival_data_AMP_Cox_AURKA_MYCN_all$High_risk <- factor(survival_data_AMP_Cox_AURKA_MYCN_all$High_risk,levels = c('Low risk','High risk','Not defined'))
survival_data_AMP_Cox_AURKA_MYCN_all$Gender<- factor(survival_data_AMP_Cox_AURKA_MYCN_all$Gender, levels = c('Not defined','Male','Female'))
survival_data_AMP_Cox_AURKA_MYCN_all$Gender<- factor(survival_data_AMP_Cox_AURKA_MYCN_all$Gender, levels = c('Male','Female','Not defined'))
survival_data_AMP_Cox_AURKA_MYCN_all$AURKA<- factor(survival_data_AMP_Cox_AURKA_MYCN_all$AURKA, levels = c('Low expression','High expression'))
survival_data_AMP_Cox_AURKA_MYCN_all$MYC1<- factor(survival_data_AMP_Cox_AURKA_MYCN_all$MYC1, levels = c('Low expression','High expression'))
survival_data_AMP_Cox_AURKA_MYCN_all$MYC2<- factor(survival_data_AMP_Cox_AURKA_MYCN_all$MYC2, levels = c('Low expression','High expression'))
survival_data_AMP_Cox_AURKA_MYCN_all$MYCN_score<- factor(survival_data_AMP_Cox_AURKA_MYCN_all$MYCN_score, levels = c('Low expression','High expression'))
survival_data_AMP_Cox_AURKA_MYCN_all$MYCN_score_PNAS<- factor(survival_data_AMP_Cox_AURKA_MYCN_all$MYCN_score_PNAS, levels = c('Low expression','High expression'))


Gcox<-coxph(Basurv~Subgroup+Age+INSS_Stage+Gender+High_risk+AURKA+MYCN_score_PNAS,data = survival_data_AMP_Cox_AURKA_MYCN_all)
Gcox<-coxph(Basurv~AURKA+Age+INSS_Stage+Gender,data = survival_data_AMP_Cox_AURKA_MYCN_all)
Gcox<-coxph(Basurv~Subgroup+High_risk,data = survival_data_AMP_Cox_AURKA_MYCN_all)

Gcox<-coxph(Basurv~MYCN_score_PNAS+Age+INSS_Stage+Gender+High_risk,data = survival_data_AMP_Cox_AURKA_MYCN_all)
Gcox<-coxph(Basurv~MYCN_score_PNAS+Age+INSS_Stage+Gender,data = survival_data_AMP_Cox_AURKA_MYCN_all)
Gcox<-coxph(Basurv~MYCN_score_PNAS+AURKA+High_risk,data = survival_data_AMP_Cox_AURKA_MYCN_all)
Gcox<-coxph(Basurv~MYCN_score_PNAS+High_risk,data = survival_data_AMP_Cox_AURKA_MYCN_all)
Gcox<-coxph(Basurv~Subgroup+High_risk,data = survival_data_AMP_Cox_AURKA_MYCN_all)
Gcox<-coxph(Basurv~AURKA+Gender,data = survival_data_AMP_Cox_AURKA_MYCN_all)
Gcox<-coxph(Basurv~INSS_Stage,data = survival_data_AMP_Cox_AURKA_MYCN_all)

colnames(survival_data_AMP_Cox_AURKA_MYCN_all)
table(survival_data_AMP_Cox_AURKA_MYCN_all$High_risk)
ggforest(Gcox,data = survival_data_AMP_Cox_AURKA_MYCN_all,cpositions = c(0.02,0.22,0.4),main = 'Hazard ratio',
         fontsize =1.8,refLabel = 'reference',noDigits = 2) #width 18 height 12, 或者7，15
?ggforest
colnames(survival_data_AMP_Cox_AURKA_MYCN_all)
install.packages('pec')
?ggforest
rownames(Clinical_test)<-Clinical_test$Sample
Clinical_test_amp_normal<-Clinical_test[c(colnames(merge_data_batch_MYCN_amp),
                                          colnames(train_1),
                                          colnames(train_2),
                                          colnames(train_3),
                                          colnames(test_1),
                                          colnames(test_2),
                                          colnames(test_3)),]
Clinical_test_amp_normal<-na.omit(Clinical_test_amp_normal)
survival_data<-merge(survival_data,Clinical_test_amp_normal,by='Sample')
for (i in 1:1530){
  if (survival_data[i,9]== 'Not defined'){
    survival_data[i,9] <-0
  } else if (survival_data[i,9]== 'No'){
    survival_data[i,9] <-1
  } else if (survival_data[i,9]== 'Yes'){
    survival_data[i,9] <-2
  }
}
table(survival_data$High.risk)
survival_data_amp<-subset(survival_data,survival_data$Group =='MYCN-AMP')
survival_data_normal<-subset(survival_data,survival_data$Group !='MYCN-AMP')
Basurv<-Surv(time=survival_data_normal$Day_Survival,event = survival_data_normal$Event)
survival_data_normal$Basurv<-with(survival_data_normal,Basurv)
Gcox<-coxph(Basurv~Group+High.risk,data = survival_data_normal)
colnames(survival_data)
class(survival_data[1,9])
ggforest(Gcox,data = survival_data_normal,cpositions = c(0.02,0.22,0.4),main = 'Hazard ratio',
         fontsize =1.8,refLabel = 'reference',noDigits = 2)

#### Prediction error curve ####
BiocManager::install('pec')
library(pec)
set.seed(123)
library(prodlim)
library(survival)
survival_data_AMP_Cox_AURKA_MYCN_all_risk<-survival_data_AMP_Cox_AURKA_MYCN_all[survival_data_AMP_Cox_AURKA_MYCN_all$High_risk != 'Not defined',]
Basurv1<-Surv(time=survival_data_AMP_Cox_AURKA_MYCN_all_risk$Day_Survival,
             event = survival_data_AMP_Cox_AURKA_MYCN_all_risk$Event)
Basurv2<-Surv(time=survival_data_AMP_Cox_AURKA_MYCN_all$Day_Survival,event = survival_data_AMP_Cox_AURKA_MYCN_all$Event)
survival_data_AMP_Cox_AURKA_MYCN_all_risk$High_risk <- factor(survival_data_AMP_Cox_AURKA_MYCN_all_risk$High_risk,levels = c('Low risk','High risk'))
Models <- list("Subgroup"=coxph(Basurv~Subgroup,x=TRUE,y=TRUE,
                                data = survival_data_AMP_Cox_AURKA_MYCN_all),
               "Age"=coxph(Basurv~Age,x=TRUE,y=TRUE,
                           data = survival_data_AMP_Cox_AURKA_MYCN_all),
               "Gender"=coxph(Basurv~Gender,x=TRUE,y=TRUE,
                              data = survival_data_AMP_Cox_AURKA_MYCN_all),
               "INSS_Stage"=coxph(Basurv~INSS_Stage,x=TRUE,y=TRUE,
                                  data = survival_data_AMP_Cox_AURKA_MYCN_all),
               "Risk"=coxph(Basurv~High_risk,x=TRUE,y=TRUE,
                            data = survival_data_AMP_Cox_AURKA_MYCN_all),
               "AURKA"=coxph(Basurv~AURKA,x=TRUE,y=TRUE,
                             data = survival_data_AMP_Cox_AURKA_MYCN_all),
               "MYCN_Score"=coxph(Basurv~MYCN_score_PNAS,x=TRUE,y=TRUE,
                                  data = survival_data_AMP_Cox_AURKA_MYCN_all)
)

Models <- list(
               "INSS_Stage"=coxph(Basurv1~INSS_Stage,x=TRUE,y=TRUE,
                                  data = survival_data_AMP_Cox_AURKA_MYCN_all_risk)
)
dim(Basurv)
dim(with(survival_data_AMP_Cox_AURKA_MYCN_all, table(INSS_Stage, Basurv)))

table(survival_data_AMP_Cox_AURKA_MYCN_all$Event)
class(survival_data_AMP_Cox_AURKA_MYCN_all_risk$INSS_Stage)

Models <- list(
  "INSS_Stage"=coxph(Basurv1~INSS_Stage,x=TRUE,y=TRUE,
                     data = survival_data_AMP_Cox_AURKA_MYCN_all_risk)
)

Models <- list("Subgroup"=coxph(Basurv1~Subgroup,x=TRUE,y=TRUE,
                                data = survival_data_AMP_Cox_AURKA_MYCN_all_risk),
               #"Age"=coxph(Basurv~Age,x=TRUE,y=TRUE,
                #           data = survival_data_AMP_Cox_AURKA_MYCN_all),
               #"Gender"=coxph(Basurv~Gender,x=TRUE,y=TRUE,
              #                data = survival_data_AMP_Cox_AURKA_MYCN_all),
               "INSS_Stage"=coxph(Basurv1~INSS_Stage,x=TRUE,y=TRUE,
                                  data = survival_data_AMP_Cox_AURKA_MYCN_all_risk),
               "Risk"=coxph(Basurv1~High_risk,x=TRUE,y=TRUE,
                            data = survival_data_AMP_Cox_AURKA_MYCN_all_risk)
              # "Subgroup+Age"=coxph(Basurv~Subgroup+Age,x=TRUE,y=TRUE,
               #                            data = survival_data_AMP_Cox_AURKA_MYCN_all),
               #"Subgroup+Gender"=coxph(Basurv~Subgroup+Gender,x=TRUE,y=TRUE,
                #                           data = survival_data_AMP_Cox_AURKA_MYCN_all),
)

PredError <- pec(object=Models,
                 formula=Surv(Day_Survival,Event)~Subgroup+Age+INSS_Stage+High_risk,
                 data=survival_data_AMP_Cox_AURKA_MYCN_all_risk,
                 exact=TRUE,
                 cens.model="cox",
                 splitMethod="none",
                 B=0,confLevel = 0.95,
                 verbose=TRUE)
plot(PredError)
plot(PredError, xlim=c(0,5000),bty='l',las=1,
     ylim = c(0, 0.3), xlab = "Days", cex.lab=3,cex.axis=3,lwd=2, smooth = T,
     legend=F)
legend('topleft',c('Reference','Subgroup','Age','Gender','INSS_Stage',
                   'High_risk'),lwd = 2,bty='n',
       col = c(1,2, 3,4,5,6),cex=1,y.intersp=0.3)
legend('topleft',c('Reference','Subgroup','INSS_Stage',
                   'High_risk'),lwd = 2,bty='n',
       col = c(1,2, 3,4),cex=1,y.intersp=0.3)

legend('topleft',c('Reference','Secreted_score','BAL_score','HIF_score','Gap_stage',
                   'Secreted_score+Gap_stage','BAL_score+Gap_stage','HIF_score+GAP_stage'),lwd = 2,bty='n',
       col = c(1,2, 3,4,5,6,7,8),cex=1.5,
       text.width=textwidths,y.intersp=0.25)  ##x/y.intersp 用来调整标题之间的间距



#### subtype marker使用meta-analysis ####

#### GSE78061 cell line ####
setwd('./GSE78061_RAW/')

library(affy)
BiocManager::install('RSQLite')

library(oligo)
library(RSQLite)
library(AnnotationDbi)
GSE78061<-ReadAffy()
GSE78061<-affy::rma(GSE78061)
GSE78061<-exprs(GSE78061)
boxplot(GSE78061)



#### CibersortX ####
CibersortxTest<-read.csv('CibersortX-Test.csv',header = T,sep = ',')
CibersortxTest<-read.xlsx('/Users/selina.weng/Downloads/6-CibersortX-Test.xlsx',sheet=4)

dim(CibersortxTest)
CibersortxTest_sub1<-CibersortxTest[CibersortxTest$Group=='Subgroup1',]
CibersortxTest_sub2<-CibersortxTest[CibersortxTest$Group=='Subgroup2',]
CibersortxTest_sub3<-CibersortxTest[CibersortxTest$Group=='Subgroup3',]


for (i in 2:9){
  print(colnames(CibersortxTest_sub1)[i])
  Covid19_1_healthy<-as.numeric(CibersortxTest_sub2[,i])
  Covid19_1_disease<-as.numeric(CibersortxTest_sub3[,i])
  
  a<-shapiro.test(Covid19_1_healthy)
  b<-shapiro.test(Covid19_1_disease)
  if (min(a$p.value,b$p.value)<0.05){
    c<-wilcox.test(Covid19_1_healthy,Covid19_1_disease) 
    print('u-test')
    #print(c$p.value)
    print(paste0('95CI:Null'))
    print(paste0('p-value: ',c$p.value))
  } else{
    print('t-test')
    d<-var.test(Covid19_1_healthy,Covid19_1_disease)
    if (d$p.value < 0.05){
      c<-t.test(Covid19_1_healthy,Covid19_1_disease,var.equal =F)
      #print(c$p.value)
    } else {
      c<-t.test(Covid19_1_healthy,Covid19_1_disease,var.equal =T)
      #print(c$p.value)
    }
    b<-c$conf.int
    print(paste0('95CI: ',round(b[1],3),'-',round(b[2],3)))
    print(paste0('p-value: ',c$p.value))
  }
}





t.test(CibersortxTest_sub1$NE_cells,CibersortxTest_sub2$NE_cells)

8*650=5200
8*513=4104
8*519=4152
8*688=5504
CibersortxTestplot<-as.data.frame(matrix(data=NA,ncol = 3,nrow=5504))
colnames(CibersortxTestplot)<-c('Celltype','Percentage',"Group")

for (i in 1:688){
  CibersortxTestplot[(8*(i-1)+1):(8*i),1]<-colnames(CibersortxTest)[2:9]
  CibersortxTestplot[(8*(i-1)+1):(8*i),2]<-as.numeric(CibersortxTest[i,2:9])
  CibersortxTestplot[(8*(i-1)+1):(8*i),3]<-CibersortxTest[i,13]
}

library(ggpubr)
my_comparisons_train <- list( c("Subgroup1", "Subgroup2"), 
                              c("Subgroup2", "Subgroup3"),
                              c("Subgroup3", "Subgroup1"))
table(CibersortxTestplot$Group)
pdf('violnnew.pdf',width =18 ,height = 12)
ggplot(CibersortxTestplot, aes(x=Celltype, y=Percentage, fill=Group))+
  geom_boxplot()+
  scale_fill_manual(values=c("#7CAE00", "#00BFC4", "#C77CFF"))+
  stat_compare_means(comparisons = my_comparisons_train,size=10)+
  theme_classic()+
  theme(
        axis.text.x = element_text( face="bold", size=30,color = 'black',angle=45,hjust = 1),
        axis.text.y = element_text( face="bold", size=30,color = 'black'),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(color="black", size=40, face="bold"),
        plot.margin = unit(c(1,1,1,1),"cm"))

dev.off()
colnames(CibersortxTestplot)
dim(train_1)
dim(train_2)
dim(train_3)



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



#### IHC survival ####

getwd()
setwd('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/')
IHC_survival<-read.xlsx('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/survival_ihc.xlsx')
fit<-survfit(Surv(Time_AURORA,Event_AURORA) ~ Group_AURORA, data=IHC_survival)

IHC_survival$Group_AURORA<-factor(IHC_survival$Group_AURORA,levels = c('High','Low'))
KMsurvival_plot<-ggsurvplot(fit,data=IHC_survival,pval = TRUE, #show p-value of log-rank test，显示log-rank分析P
                            #conf.int = TRUE, #添加置信区间
                            pval.size=10,
                            legend.labs =  c("AURORA-low",'AURORA-high'),
                            legend.title='', 
                            xlab = "Time in Day",   ###  customize X axis label.自定义x的标time in years
                            #xlim=c(0,50),
                            break.x.by=10, ###改变坐标轴间距
                            ylab=paste0('Overall survival Days'),
                            #ylab=paste0('Overall peak day'),
                            #ylab=paste0('SRAS−CoV−2 RNA +'),
                            surv.median.line = "hv", #添加中位生存时间的线
                            palette = c( "blue","red"), ###  自定义颜色
                            #font.main = c(16, "bold", "darkblue"),
                            font.x = 40, # X 轴
                            font.y = 35, # c(14, "bold.italic", "darkred"), y 轴特征
                            font.tickslab = 30,# c(12, "plain", "darkgreen"), 坐标轴大小
                            #conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样子
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)的形式展示risk table
                            risk.table.y.text.col = T,###  colour risk table text annotations.
                            risk.table.y.text = F,###  show bars instead of names in text annotations in legend of risk table.不显示legend名字
                            tables.y.text=F,
                            tables.x.text =F,
                            #risk.table.title="My title", ## 改变 title
                            fontsize=5, ## 表格中数据大小
                            ncensor.plot = F, #我这里不显TRUE就显示
                            #tables.theme=theme_cleantable(), # 取消table边框
                            ggtheme = theme_classic()#绘图主题
)

KMsurvival_plot<-ggsurvplot(fit,data=IHC_survival,pval = TRUE, #show p-value of log-rank test，显示log-rank分析得到的P值
                            #conf.int = TRUE, #添加置信区间
                            pval.size=10,
                            legend.labs =  c('MYCN-negative',"MYCN-positive"),
                            legend.title='', 
                            xlab = "Time in Day",   ###  customize X axis label.自定义x的标time in years
                            #xlim=c(0,50),
                            break.x.by=10, ###改变坐标轴间距
                            ylab=paste0('Overall survival Days'),
                            #ylab=paste0('Overall peak day'),
                            #ylab=paste0('SRAS−CoV−2 RNA +'),
                            surv.median.line = "hv", #添加中位生存时间的线
                            palette = c( "blue","red"), ###  自定义颜色
                            #font.main = c(16, "bold", "darkblue"),
                            font.x = 40, # X 轴
                            font.y = 35, # c(14, "bold.italic", "darkred"), y 轴特征
                            font.tickslab = 30,# c(12, "plain", "darkgreen"), 坐标
                            #conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样子
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)的形式展示risk table
                            risk.table.y.text.col = T,###  colour risk table text annotations.
                            risk.table.y.text = F,###  show bars instead of names in text annotations in legend of risk table.不显示legend名字
                            tables.y.text=F,
                            tables.x.text =F,
                            #risk.table.title="My title", ## 改变 title名字
                            fontsize=5, ## 表格中数据大小
                            ncensor.plot = F, #我这里不显示删失的TRUE就显示
                            #tables.theme=theme_cleantable(), # 取消table
                            ggtheme = theme_classic()#绘图主题
)
KMsurvival_plot$plot<-KMsurvival_plot$plot+
  theme(legend.text = element_text(size = 30),plot.margin = unit(c(2,2,0,2), "cm"))

KMsurvival_plot$table<-KMsurvival_plot$table+labs(x = NULL, y = NULL)+theme(axis.text=element_text(size=30),
                                                                            axis.title=element_text(size=30),
                                                                            
                                                                            legend.text = element_text(size = 25),
                                                                            plot.title = element_text(size=25),
                                                                            plot.margin = unit(c(0,2,2,2), "cm")) # 改变numer at risk 大小)

ggsave(file=paste0('IHC_AURORA_survival.pdf'),width = 15,height = 12,print(KMsurvival_plot),
       onefile=FALSE)
table(IHC_survival$Group_AURORA)
getwd()

urvival_data[colnames(train_1),4]<-'Subgroup1' ## 6是根据Group所在的列
survival_data[colnames(train_2),4]<-'Subgroup2' ## 6是根据Group所在的列
survival_data[colnames(train_3),4]<-'Subgroup3' ## 6是根据Group所在的列
survival_data_AMP_Cox<-merge(survival_data_AMP,Clinical_test_normal,by='Sample')
sub2_up_pos_removeMYCN_matrix_t$Sample<-rownames(sub2_up_pos_removeMYCN_matrix_t)
survival_data_AMP_Cox<-merge(survival_data_AMP_Cox,sub2_up_pos_removeMYCN_matrix_t,by='Sample')
survival_data_AMP_Cox_surv<-survival_data_AMP_Cox[,c(1,2,3,10:125)]
IHC_survival<-read.csv('/Users/yihuawang/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/AURKA/Survival_median_AURKA_train+testall+AMP.csv')
Basurv<-Surv(time=IHC_survival$Day_Survival,event = IHC_survival$Event)
Basurv<-Surv(time=IHC_survival$Time_AURORA,event = IHC_survival$Event_AURORA)

IHC_survival$Basurv<-with(IHC_survival,Basurv)
Gcox<-coxph(Basurv~HIF1A_1,data=IHC_survival)
Gsum<-summary(Gcox)

round(Gsum$coefficients[,2],3)
Gsum$coefficients[,5]
paste0(round(Gsum$conf.int[,3:4],2),collapse = '-')


UniCox<-function(x){
  FML<-as.formula(paste0('Basurv~',x))
  Gcox<-coxph(Basurv~Group_MYCN,data=IHC_survival)
  Gsum<-summary(Gcox)
  
  HR<-round(Gsum$coefficients[,2],3)
  Pvalue<-Gsum$coefficients[,5]
  CI<-paste0(round(Gsum$conf.int[,3:4],2),collapse = '-')
  Unicox<-data.frame('Characteristics'=x,
                     'Hazard Ratio'=HR,
                     'CI95'=CI,
                     'P value'=Pvalue)
  return(Unicox)
}

#### ROC for svm model ####
predicted.classes_lasso<-as.data.frame(predicted.classes_lasso)
colnames(predicted.classes_lasso)<-c(1,2,3)
multiclass.roc(response = test_GM_zscore$Group, predictor =as.data.frame(predicted.classes_lasso),
               plot=T,print.auc=TRUE)

predicted.classes_lasso<-as.data.frame(predicted.classes_lasso)
colnames(predicted.classes_lasso)<-c(1,2,3)
multiclass.roc(response = test_GM_zscore$Group, predictor =as.data.frame(predicted.classes_lasso),
               plot=T,print.auc=TRUE)
30*2+18*2=96
24/96*18*2
24/96*30
167*2+39*2+98+37
# Test
max(attr(pred_zscore_removesomegene_score, "probabilities")[1:5,])
pred_zscore_removesomegene[1:5,]
library(pROC)
pdf('Test_AUC.pdf',width = 5,height = 5)
par(mar=c(4,4,4,4))
multiclass.roc(response = test_GM_zscore$Group, predictor =attr(pred_zscore_removesomegene_score, "probabilities"),
                               plot=T,print.auc=TRUE)
dev.off()
?multiclass.roc
#ggroc(roc_svm_test$rocs)
#rs <- roc_svm_test[['rocs']]

# Fudan 
pdf('FUDAN_AUC.pdf',width = 5,height = 5)
par(mar=c(4,4,4,4))
multiclass.roc(response = train_test_GM_GSe49710$Group, predictor =attr(pred_zscore_core_score, "probabilities"),
               plot=T,print.auc=TRUE)
dev.off()

# TARGET microarray
pdf('TARGET_microarray_AUC.pdf',width = 5,height = 5)
par(mar=c(4,4,4,4))
multiclass.roc(response = Target_clinical$Subgroup, predictor =attr(pred_Target_micro_score, "probabilities"),
               plot=T,print.auc=TRUE)
dev.off()
table(Target_clinical_RNA_MYCNnormal$Subgroup)
# TARGET RNA-seq from silhou.R
Target_clinical_RNA_MYCNnormal

pdf('TARGET_RNA_AUC.pdf',width = 5,height = 5)
par(mar=c(4,4,4,4))
multiclass_roc<-multiclass.roc(response = Target_clinical_RNA_MYCNnormal$Subgroup, predictor =attr(pred_Target_score , "probabilities"),
               plot=T,print.auc=TRUE)
dev.off()
table(pred_Target_zscore$pred_Target_zscore,Target_clinical_RNA_MYCNnormal$Subgroup)
pred_Target_zscore
ci(multiclass_roc)

save(predicted.classes_lasso,pred_zscore_removesomegene_score,pred_zscore_core_score,
     pred_Target_micro_score,pred_Target_score,
     file='ML_predicted_probabilities.Rdata')
dim(attr(pred_zscore_core_score, "probabilities"))


library(ggplot2)
library(ggExtra)
BiocManager::install('ggExtra')
Probability_Test<-as.data.frame(attr(pred_zscore_removesomegene_score, "probabilities")) # Test
Probability_Test$Sample<-rownames(test_GM)
Probability_Test$Group<-pred_zscore_removesomegene
Probability_Test$Cluster<-test_GM_zscore$Group

Probability_Test<-as.data.frame(attr(pred_Target_micro_score, "probabilities")) # TARGET array
Probability_Test$Sample<-Target_clinical$TARGET.USI
Probability_Test$Group<-pred_Target_micro_zscore$pred_Target_micro_zscore
Probability_Test$Cluster<-Target_clinical$Subgroup

Probability_Test<-as.data.frame(attr(pred_Target_score, "probabilities")) # TARGET RNA
Probability_Test$Sample<-Target_clinical_RNA_MYCNnormal$Sample_2
Probability_Test$Group<-pred_Target_zscore$pred_Target_zscore
Probability_Test$Cluster<-Target_clinical_RNA_MYCNnormal$Subgroup

Probability_Test<-as.data.frame(attr(pred_zscore_core_score, "probabilities")) # GSE49711 RNA
Probability_Test$Sample<-rownames(train_test_GM_GSe49710)
Probability_Test$Group<-pred_zscore_core
Probability_Test$Cluster<-train_test_GM_GSe49710$Group


colnames(Probability_Test)[1:3]<-c('Subgroup1','Subgroup2','Subgroup3')
write.table(Probability_Test,
            'Probability_Test.txt',sep = '\t')
write.table(Probability_Test,
            'Probability_TARGET_Microarray.txt',sep = '\t')
write.table(Probability_Test,
            'Probability_TARGET_RNAseq.txt',sep = '\t')
write.table(Probability_Test,
            'Probability_GSE49711_RNAseq.txt',sep = '\t')

p<-ggplot(Probability_Test, aes(x=Subgroup2, y=Subgroup3, color=Group)) + 
  geom_point(size=3)+
  theme_bw()+
  scale_color_manual(values=c('#7CAE00','#00BFC4', '#C77CFF'))+
  theme(legend.position = "none",
        axis.text.x = element_text( size=18,color = 'black',hjust = 0.5,vjust=0),
        axis.text.y = element_text( size=18,color = 'black'),
        #axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(color="black", size=28),
        axis.title.x = element_text(color="black", size=28),
        plot.margin = unit(c(4,4,4,4), "cm"))
p1 <- ggMarginal(p, type="densigram", groupColour = T,groupFill=T,
                 size = 5)
p1
ggsave(file=paste0('Test_Probability2.pdf'),width = 6,height = 6,print(p1),
       onefile=FALSE)


#####
Fudan_sankeywheel<-read.xlsx('GEO+Array/Consensus/Kmean/Figure folder/Figure8 Prediction/Fudan_Cluster_Prediction.xlsx',
                                           sheet = 1)
dev.new()
sankeywheel(from = Fudan_sankeywheel$from,
            to = Fudan_sankeywheel$to,
            weight =Fudan_sankeywheel$weight,
            type = "sankey", 
            width = "100%",
            title = "")
#### Comparison with public subgroups ####

Zhou_subgroup<-as.data.frame(matrix(data =NA,ncol=2,nrow=3085))
colnames(Zhou_subgroup)<-c('Zhou_SampleID','Zhou_Subgroup')
Zhou_subgroup$Zhou_SampleID[1:1034]<-colnames(train_1)
Zhou_subgroup$Zhou_Subgroup[1:1034]<-'Subgroup1'

Zhou_subgroup$Zhou_SampleID[1035:1593]<-colnames(train_2)
Zhou_subgroup$Zhou_Subgroup[1035:1593]<-'Subgroup2'
Zhou_subgroup$Zhou_SampleID[1594:2160]<-colnames(train_3)
Zhou_subgroup$Zhou_Subgroup[1594:2160]<-'Subgroup3'
Zhou_subgroup$Zhou_SampleID[2161:2538]<-colnames(test_1)
Zhou_subgroup$Zhou_Subgroup[2161:2538]<-'Subgroup1'
Zhou_subgroup$Zhou_SampleID[2539:2806]<-colnames(test_2)
Zhou_subgroup$Zhou_Subgroup[2539:2806]<-'Subgroup2'
Zhou_subgroup$Zhou_SampleID[2807:3085]<-colnames(test_3)
Zhou_subgroup$Zhou_Subgroup[2807:3085]<-'Subgroup3'
rownames(Zhou_subgroup)<-Zhou_subgroup$Zhou_SampleID
write.table(Zhou_subgroup,'Zhou_subgroup_infromation.txt',sep = '\t')

# GSE85047 Andrea Califano Cancer Discov #
library(openxlsx)
getwd()
GSE85047_Andrea<-read.xlsx('GEO+Array/Consensus/Kmean/Comparison_Public_Subgroup/Andrea_Califano_Cancer_Discovery_TARGET_GSE85047.xlsx',
                           sheet = 2)
GSE85047_Zhou<-read.xlsx('GEO+Array/Consensus/Kmean/Comparison_Public_Subgroup/Andrea_Califano_Cancer_Discovery_TARGET_GSE85047.xlsx',
                         sheet = 3)
GSE85047_Zhou<-merge(GSE85047_Zhou,GSE85047_Andrea,by='Andrea_SampleID')
GSE85047_Zhou<-na.omit(GSE85047_Zhou)
GSE85047_Zhou<-GSE85047_Zhou[GSE85047_Zhou$`MYCN-status`==0,]
GSE85047_Zhou<-merge(GSE85047_Zhou,Zhou_subgroup,by='Zhou_SampleID')
table(GSE85047_Zhou$Andrea_Subgroup,GSE85047_Zhou$Zhou_Subgroup)
write.csv(GSE85047_Zhou,'GSE85047_Zhou.csv')

### TARGET Andrea Califano Cancer Discov ###
TARGET_Andrea<-read.xlsx('GEO+Array/Consensus/Kmean/Comparison_Public_Subgroup/Andrea_Califano_Cancer_Discovery_TARGET_GSE85047.xlsx',
                         sheet = 1)
TARGET_NBL_SampleMatrix<-read.xlsx('GEO+Array/Consensus/Kmean/Comparison_Public_Subgroup/TARGET_NBL_SampleMatrix_Discovery_20180111.xlsx',
                                   sheet = 1)
TARGET_Andrea<-merge(TARGET_Andrea,TARGET_NBL_SampleMatrix,by='Sample.ID')
TARGET_Andrea<-TARGET_Andrea[TARGET_Andrea$`MYCN-status`=='0',]
TARGET_Andrea<-TARGET_Andrea[,1:6]

colnames(TARGET_Andrea)[6]<-'TARGET.USI'

TARGET_Andrea<-merge(TARGET_Andrea,Target_clinical,by='TARGET.USI') # Target_Clinical from silhou.R

write.table(TARGET_Andrea,'TARGET_Andrea.txt',sep = '\t')
table(TARGET_Andrea$SUBTYPE,TARGET_Andrea$Subgroup)


library(sankeywheel)
getwd()
GSE85047_Andrea_sankeywheel<-read.xlsx('GEO+Array/Consensus/Kmean/Figure folder/Comparison_Public_Subgroup/Andrea_Califano_Cancer_Discovery_TARGET_GSE85047.xlsx',
                                       sheet = 6)
GSE85047_Andrea_sankeywheel<-read.xlsx('/Volumes/My Passport/2021/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/Figure folder/Comparison_Public_Subgroup/Andrea_Califano_Cancer_Discovery_TARGET_GSE85047.xlsx',
                                       sheet = 5)
GSE49710_Westermann<-read.xlsx('GEO+Array/Consensus/Kmean/Figure folder/Comparison_Public_Subgroup/Westermann_Nature_Cancer_GSE49710.xlsx',
                     sheet = 3)

GSE49710_Westermann<-read.xlsx('/Volumes/My Passport/2021/Neuroblastoma/Hu/GEO+Array/Consensus/Kmean/Figure folder/Comparison_Public_Subgroup/Westermann_Nature_Cancer_GSE49710.xlsx',
                               sheet = 3)


GSE49710_Westermann<-merge(GSE49710_Westermann,Zhou_subgroup,by='Zhou_SampleID')
write.csv(GSE49710_Westermann,'GSE49710_Westermann.csv')
table(GSE49710_Westermann$Signature,GSE49710_Westermann$Zhou_Subgroup)
GSE49710_Westermann_sankeywheel<-read.xlsx('GEO+Array/Consensus/Kmean/Comparison_Public_Subgroup/Westermann_Nature_Cancer_GSE49710.xlsx',
                                           sheet = 4)
dev.new()
sankeywheel(from = GSE85047_Andrea_sankeywheel$fro,
            to = GSE85047_Andrea_sankeywheel$to,
            weight =GSE85047_Andrea_sankeywheel$weight,
            type = "sankey", 
            width = "100%",
            title = "")
#### 

{GSE49710_Westermann_mes_sub3<-GSE49710_Westermann[GSE49710_Westermann$Signature=='MES',]
GSE49710_Westermann_mes_sub3<-GSE49710_Westermann_mes_sub3[GSE49710_Westermann_mes_sub3$Zhou_Subgroup=='Subgroup3',]

GSE49710_Westermann_HR_sub1<-GSE49710_Westermann[GSE49710_Westermann$Signature=='MNA-HR',]
GSE49710_Westermann_HR_sub1<-GSE49710_Westermann_HR_sub1[GSE49710_Westermann_HR_sub1$Zhou_Subgroup=='Subgroup1',]


GSE49710_Westermann_HR_sub2<-GSE49710_Westermann[GSE49710_Westermann$Signature=='MNA-HR',]
GSE49710_Westermann_HR_sub2<-GSE49710_Westermann_HR_sub2[GSE49710_Westermann_HR_sub2$Zhou_Subgroup=='Subgroup2',]
GSE49710_Westermann_mes_HR<-rbind(GSE49710_Westermann_mes_sub3,
                                  GSE49710_Westermann_HR_sub1,
                                  GSE49710_Westermann_HR_sub2)
}
#
{
GSE49710_Westermann_LR_sub1<-GSE49710_Westermann[GSE49710_Westermann$Signature=='MNA-LR',]
GSE49710_Westermann_HR_sub1<-GSE49710_Westermann_LR_sub1[GSE49710_Westermann_LR_sub1$Zhou_Subgroup=='Subgroup1',]
GSE49710_Westermann_LR_sub2<-GSE49710_Westermann[GSE49710_Westermann$Signature=='MNA-LR',]
GSE49710_Westermann_LR_sub2<-GSE49710_Westermann_LR_sub1[GSE49710_Westermann_LR_sub1$Zhou_Subgroup=='Subgroup2',]
GSE49710_Westermann_LR_sub3<-GSE49710_Westermann[GSE49710_Westermann$Signature=='MNA-LR',]
GSE49710_Westermann_LR_sub3<-GSE49710_Westermann_LR_sub3[GSE49710_Westermann_LR_sub3$Zhou_Subgroup=='Subgroup3',]
GSE49710_Westermann_LR<-rbind(
                                  GSE49710_Westermann_LR_sub1,
                                  GSE49710_Westermann_LR_sub2,
                                  GSE49710_Westermann_LR_sub3)
}
#
{
  GSE49710_Sub3_LR<-GSE49710_Westermann[GSE49710_Westermann$Zhou_Subgroup=='Subgroup3',]
  GSE49710_Sub3_LR<-GSE49710_Sub3_LR[GSE49710_Sub3_LR$Signature=='MNA-LR',]
  GSE49710_Sub3_HR<-GSE49710_Westermann[GSE49710_Westermann$Zhou_Subgroup=='Subgroup3',]
  GSE49710_Sub3_HR<-GSE49710_Sub3_HR[GSE49710_Sub3_HR$Signature=='MNA-HR',]
  GSE49710_Sub3_MES<-GSE49710_Westermann[GSE49710_Westermann$Zhou_Subgroup=='Subgroup3',]
  GSE49710_Sub3_MES<-GSE49710_Sub3_MES[GSE49710_Sub3_MES$Signature=='MES',]
  GSE49710_Sub3_LRHRMES<-rbind(
    GSE49710_Sub3_LR,
    GSE49710_Sub3_HR,
    GSE49710_Sub3_MES)
  GSE49710_Sub3_LRHRMES$Signature<-factor(GSE49710_Sub3_LRHRMES$Signature,levels = c("MNA-LR",'MNA-HR','MES'))
}
#Sub1--LR/HR/MES
{
  GSE49710_Sub1_LR<-GSE49710_Westermann[GSE49710_Westermann$Zhou_Subgroup=='Subgroup1',]
  GSE49710_Sub1_LR<-GSE49710_Sub1_LR[GSE49710_Sub1_LR$Signature=='MNA-LR',]
  GSE49710_Sub1_HR<-GSE49710_Westermann[GSE49710_Westermann$Zhou_Subgroup=='Subgroup1',]
  GSE49710_Sub1_HR<-GSE49710_Sub1_HR[GSE49710_Sub1_HR$Signature=='MNA-HR',]
  GSE49710_Sub1_MES<-GSE49710_Westermann[GSE49710_Westermann$Zhou_Subgroup=='Subgroup1',]
  GSE49710_Sub1_MES<-GSE49710_Sub1_MES[GSE49710_Sub1_MES$Signature=='MES',]
  GSE49710_Sub1_LRHRMES<-rbind(
    GSE49710_Sub1_LR,
    GSE49710_Sub1_HR,
    GSE49710_Sub1_MES)
  GSE49710_Sub1_LRHRMES$Signature<-factor(GSE49710_Sub1_LRHRMES$Signature,levels = c("MNA-LR",'MNA-HR','MES'))
}

{
  GSE49710_Sub1_LR<-GSE49710_Westermann[GSE49710_Westermann$Zhou_Subgroup=='Subgroup2',]
  GSE49710_Sub1_LR<-GSE49710_Sub1_LR[GSE49710_Sub1_LR$Signature=='MNA-LR',]
  GSE49710_Sub1_HR<-GSE49710_Westermann[GSE49710_Westermann$Zhou_Subgroup=='Subgroup2',]
  GSE49710_Sub1_HR<-GSE49710_Sub1_HR[GSE49710_Sub1_HR$Signature=='MNA-HR',]
  GSE49710_Sub1_MES<-GSE49710_Westermann[GSE49710_Westermann$Zhou_Subgroup=='Subgroup2',]
  GSE49710_Sub1_MES<-GSE49710_Sub1_MES[GSE49710_Sub1_MES$Signature=='MES',]
  GSE49710_Sub1_LRHRMES<-rbind(
    GSE49710_Sub1_LR,
    GSE49710_Sub1_HR,
    GSE49710_Sub1_MES)
  GSE49710_Sub1_LRHRMES$Signature<-factor(GSE49710_Sub1_LRHRMES$Signature,levels = c("MNA-LR",'MNA-HR','MES'))
}
fit<-survfit(Surv(Day_survival,Event) ~ Zhou_Subgroup, data=GSE49710_Westermann_mes_HR)
fit<-survfit(Surv(Day_survival,Event) ~ Zhou_Subgroup, data=GSE49710_Westermann_LR)
fit<-survfit(Surv(Day_survival,Event) ~ Signature, data=GSE49710_Sub3_LRHRMES)
fit<-survfit(Surv(Day_survival,Event) ~ Signature, data=GSE49710_Sub1_LRHRMES)
factor(GSE49710_Sub1_LRHRMES$Signature)

Hr15,lr24,mes63
KMsurvival_plot<-ggsurvplot(fit,data=GSE49710_Sub1_LRHRMES,pval = TRUE, #show p-value of log-rank test，显示log-rank分析得到的P值
                            #conf.int = TRUE, #添加置信区间
                            pval.size=10,
                            #legend.labs =  c("Subgroup1",'Subgroup2','Subgroup3'),
                            legend.labs =  c("LR",'HR','MES'),
                            
                            legend.title='', 
                            xlab = "Time in Day",   ###  customize X axis label.自定义x的标time in years
                            #xlim=c(0,50),
                            break.x.by=1000, ###改变坐标轴间距
                            ylab=paste0('Overall survival Days'),
                            #ylab=paste0('Overall peak day'),
                            #ylab=paste0('SRAS−CoV−2 RNA +'),
                            surv.median.line = "hv", #添加中位生存时间的线
                            palette = c('#7CAE00','#00BFC4', '#C77CFF'), ###  自定义颜色
                            #font.main = c(16, "bold", "darkblue"),
                            font.x = 40, # X 轴
                            font.y = 35, # c(14, "bold.italic", "darkred"), y 轴特征
                            font.tickslab = 30,# c(12, "plain", "darkgreen"), 坐标轴大小
                            #conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样子
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)的形式展示risk table
                            risk.table.y.text.col = T,###  colour risk table text annotations.
                            risk.table.y.text = F,###  show bars instead of names in text annotations in legend of risk table.不显示legend名字
                            tables.y.text=F,
                            tables.x.text =F,
                            #risk.table.title="My title", ## 改变 title名字
                            fontsize=5, ## 表格中数据大小
                            ncensor.plot = F, #我这里不显示删失的TRUE就显示
                            #tables.theme=theme_cleantable(), # 取消table边框
                            ggtheme = theme_classic()#绘图主题
)
KMsurvival_plot$plot<-KMsurvival_plot$plot+
  theme(legend.text = element_text(size = 30),plot.margin = unit(c(2,2,0,2), "cm"))

KMsurvival_plot$table<-KMsurvival_plot$table+labs(x = NULL, y = NULL)+theme(axis.text=element_text(size=30),
                                                                            axis.title=element_text(size=30),
                                                                            
                                                                            legend.text = element_text(size = 25),
                                                                            plot.title = element_text(size=25),
                                                                            plot.margin = unit(c(0,2,2,2), "cm")) # 改变numer at risk 大小)

ggsave(file=paste0('GSE49710_Sub1_LRHRMES.pdf'),width = 15,height = 12,print(KMsurvival_plot),
       onefile=FALSE)

# GSE85047 Andrea Califano Cancer Discov #

GSE85047_Andrea_11q_sub1<-GSE85047_Andrea_sankeywheel[GSE85047_Andrea_sankeywheel$Andrea_Subgroup=='11q-LOH',]
GSE85047_Andrea_11q_sub1<-GSE85047_Andrea_11q_sub1[GSE85047_Andrea_11q_sub1$Zhou_Subgroup=='Subgroup1',]

GSE85047_Andrea_11q_sub2<-GSE85047_Andrea_sankeywheel[GSE85047_Andrea_sankeywheel$Andrea_Subgroup=='11q-LOH',]
GSE85047_Andrea_11q_sub2<-GSE85047_Andrea_11q_sub2[GSE85047_Andrea_11q_sub2$Zhou_Subgroup=='Subgroup2',]

GSE85047_Andrea_11q_sub3<-GSE85047_Andrea_sankeywheel[GSE85047_Andrea_sankeywheel$Andrea_Subgroup=='11q-LOH',]
GSE85047_Andrea_11q_sub3<-GSE85047_Andrea_11q_sub3[GSE85047_Andrea_11q_sub3$Zhou_Subgroup=='Subgroup3',]
#MES
GSE85047_Andrea_11q_sub1<-GSE85047_Andrea_sankeywheel[GSE85047_Andrea_sankeywheel$Andrea_Subgroup=='MES',]
GSE85047_Andrea_11q_sub1<-GSE85047_Andrea_11q_sub1[GSE85047_Andrea_11q_sub1$Zhou_Subgroup=='Subgroup1',]

GSE85047_Andrea_11q_sub2<-GSE85047_Andrea_sankeywheel[GSE85047_Andrea_sankeywheel$Andrea_Subgroup=='MES',]
GSE85047_Andrea_11q_sub2<-GSE85047_Andrea_11q_sub2[GSE85047_Andrea_11q_sub2$Zhou_Subgroup=='Subgroup2',]

GSE85047_Andrea_11q_sub3<-GSE85047_Andrea_sankeywheel[GSE85047_Andrea_sankeywheel$Andrea_Subgroup=='MES',]
GSE85047_Andrea_11q_sub3<-GSE85047_Andrea_11q_sub3[GSE85047_Andrea_11q_sub3$Zhou_Subgroup=='Subgroup3',]

GSE85047_Andrea_11q<-rbind(GSE85047_Andrea_11q_sub1,
                           GSE85047_Andrea_11q_sub2,GSE85047_Andrea_11q_sub3)
colnames(GSE85047_Andrea_11q)[1]<-'Sample'
GSE85047_Andrea_11q<-merge(GSE85047_Andrea_11q,survival_data,by='Sample')

TARGET_Andrea_11q_sub1<-GSE85047_Andrea_sankeywheel[GSE85047_Andrea_sankeywheel$Andrea_SUBTYPE=='11q-LOH',]
TARGET_Andrea_11q_sub1<-TARGET_Andrea_11q_sub1[TARGET_Andrea_11q_sub1$Zhou_Subgroup=='Subgroup1',]

TARGET_Andrea_11q_sub2<-GSE85047_Andrea_sankeywheel[GSE85047_Andrea_sankeywheel$Andrea_SUBTYPE=='11q-LOH',]
TARGET_Andrea_11q_sub2<-TARGET_Andrea_11q_sub2[TARGET_Andrea_11q_sub2$Zhou_Subgroup=='Subgroup2',]

TARGET_Andrea_11q_sub3<-GSE85047_Andrea_sankeywheel[GSE85047_Andrea_sankeywheel$Andrea_SUBTYPE=='11q-LOH',]
TARGET_Andrea_11q_sub3<-TARGET_Andrea_11q_sub3[TARGET_Andrea_11q_sub3$Zhou_Subgroup=='Subgroup3',]

TARGET_Andrea_11q<-rbind(TARGET_Andrea_11q_sub1,
                         TARGET_Andrea_11q_sub2,TARGET_Andrea_11q_sub3)

fit<-survfit(Surv(Day_Survival,Event) ~ Zhou_Subgroup, data=GSE85047_Andrea_11q)
fit<-survfit(Surv(Overall.Survival.Time.in.Days,Vital.Status.1) ~ Zhou_Subgroup, data=TARGET_Andrea_11q)


KMsurvival_plot<-ggsurvplot(fit,data=TARGET_Andrea_11q,pval = TRUE, #show p-value of log-rank test，显示log-rank分析得到的P值
                            #conf.int = TRUE, #添加置信区间
                            pval.size=10,
                            legend.labs =  c("Subgroup1",'Subgroup2','Subgroup3'),
                            
                            legend.title='', 
                            xlab = "Time in Day",   ###  customize X axis label.自定义x的标time in years
                            #xlim=c(0,50),
                            break.x.by=1000, ###改变坐标轴间
                            ylab=paste0('Overall survival Days'),
                            #ylab=paste0('Overall peak day'),
                            #ylab=paste0('SRAS−CoV−2 RNA +'),
                            surv.median.line = "hv", #添加中位生存时间的线
                            palette = c('#7CAE00','#00BFC4', '#C77CFF'), ###  自定义颜色
                            #font.main = c(16, "bold", "darkblue"),
                            font.x = 40, # X 轴
                            font.y = 35, # c(14, "bold.italic", "darkred"), y 轴特征
                            font.tickslab = 30,# c(12, "plain", "darkgreen"), 坐标轴大小
                            #conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样子
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)的形式展示risk table
                            risk.table.y.text.col = T,###  colour risk table text annotations.
                            risk.table.y.text = F,###  show bars instead of names in text annotations in legend of risk table.不显示legend名字
                            tables.y.text=F,
                            tables.x.text =F,
                            #risk.table.title="My title", ## 改变 title名字
                            fontsize=5, ## 表格中数据大小
                            ncensor.plot = F, #我这里不显TRUE就显示
                            #tables.theme=theme_cleantable(), # 取消table边框
                            ggtheme = theme_classic()#绘图主题
)
KMsurvival_plot$plot<-KMsurvival_plot$plot+
  theme(legend.text = element_text(size = 30),plot.margin = unit(c(2,2,0,2), "cm"))

KMsurvival_plot$table<-KMsurvival_plot$table+labs(x = NULL, y = NULL)+theme(axis.text=element_text(size=30),
                                                                            axis.title=element_text(size=30),
                                                                            
                                                                            legend.text = element_text(size = 25),
                                                                            plot.title = element_text(size=25),
                                                                            plot.margin = unit(c(0,2,2,2), "cm")) # 改变numer at risk 大小)

ggsave(file=paste0('Test_Km_3_TARGET_Andrea_11q.pdf'),width = 15,height = 12,print(KMsurvival_plot),
       onefile=FALSE)

# Ruidong_Group Assignment #
Ruidong_GSE49710<-read.xlsx('GEO+Array/Consensus/Kmean/Comparison_Public_Subgroup/Ruidong_Cancercell_GSE49710.xlsx',
                            sheet = 1)
Ruidong_GSE49710<-merge(Ruidong_GSE49710,Zhou_subgroup,by='Zhou_SampleID')
write.csv(Ruidong_GSE49710,'Ruidong_GSE49710.csv')
table(Ruidong_GSE49710$group,Ruidong_GSE49710$Zhou_Subgroup)
Ruidong_GSE49710<-read.xlsx('GEO+Array/Consensus/Kmean/Comparison_Public_Subgroup/Ruidong_Cancercell_GSE49710.xlsx',
                            sheet = 3)
sankeywheel(from = Ruidong_GSE49710$from,
            to = Ruidong_GSE49710$to,
            weight =Ruidong_GSE49710$weight,
            type = "sankey", 
            width = "100%",
            title = "")


#### Rogier_Versteeg_NG_ADRN_MES ####

ADRN_siganture<-read.xlsx('Rogier_Versteeg_NG_ADRN_MES.xlsx')
ADRN_siganture<-ADRN_siganture[ADRN_siganture$Group=='ADRN',]
ADRN_siganture<-ADRN_siganture[ADRN_siganture$Group!='ADRN',]

list(MYCN_signature_po$Gene)
gsva_ADRN_matrix<-as.data.frame(matrix(data=NA,
                                       nrow=3084,ncol=7))

gsva_ADRN_matrix[1:714,2]<-t(gsva(as.matrix(merge_data_batch_MYCN_amp),list(ADRN_siganture$Gene) , mx.diff=1,kcdf="Gaussian",
                                 method= "ssgsea",ssgsea.norm=F))

gsva_ADRN_matrix[715:1748,2]<-t(gsva(as.matrix(train_1_sil),list(ADRN_siganture$Gene) , mx.diff=1,kcdf="Gaussian",
                                    method= "ssgsea",ssgsea.norm=F))
gsva_ADRN_matrix[1749:2126,2]<-t(gsva(as.matrix(test_1_sil),list(ADRN_siganture$Gene), mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_ADRN_matrix[2127:2442,2]<-t(gsva(as.matrix(train_2_sil),list(ADRN_siganture$Gene) , mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_ADRN_matrix[2443:2550,2]<-t(gsva(as.matrix(test_2_sil),list(ADRN_siganture$Gene), mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_ADRN_matrix[2551:2920,2]<-t(gsva(as.matrix(train_3_sil),list(ADRN_siganture$Gene), mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_ADRN_matrix[2921:3084,2]<-t(gsva(as.matrix(test_3_sil),list(ADRN_siganture$Gene), mx.diff=1,kcdf="Gaussian",
                                     method= "ssgsea",ssgsea.norm=F))
gsva_ADRN_matrix[,1]<-c(colnames(merge_data_batch_MYCN_amp),
                        colnames(train_1_sil),colnames(test_1_sil),
                        colnames(train_2_sil),colnames(test_2_sil),
                        colnames(train_3_sil),colnames(test_3_sil))
colnames(gsva_MYC_matrix)
gsva_ADRN_matrix[,3:6]<-gsva_MYC_matrix[,6:9]
colnames(gsva_ADRN_matrix)<-c('SampleID','ADRN_score','Group','MYCN_status','MYCN_subgroup','Train_Test','REST')
gsva_ADRN_matrix_test<-gsva_ADRN_matrix[gsva_ADRN_matrix$Train_Test=='Train',]

my_comparisons_train <- list( c("Subgroup1", "Subgroup2"), 
                              c("Subgroup2", "Subgroup3"),
                              c("Subgroup3", "Subgroup1"))
ggviolin(gsva_ADRN_matrix_test, x = "MYCN_subgroup", y = 'ADRN_score',
         fill = "MYCN_subgroup",
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons_train,size=10)+
  ylab("ADRN score\n") +
  xlab('')+
  #scale_y_continuous(breaks=seq(0, 1, 0.1))+
  #coord_cartesian(ylim = c(-2500,3000))+
  #geom_hline(yintercept = median(gsva_MYC_matrix[Purity_Test$Group=='MYCN_AMP',]$ImmuneScore), 
  #          size=2,linetype = 2)+
  theme_classic()+
  theme(legend.position = "none",
        axis.text.x = element_text( face="bold", size=30,color = 'black',hjust = 0.5,vjust=0),
        axis.text.y = element_text( face="bold", size=30,color = 'black'),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid"),
        axis.title.y = element_text(color="black", size=40, face="bold"),
        plot.margin = unit(c(1,1,1,1),"cm"))
ggsave('Train_ADRN_score.pdf',width = 12,height = 12) 

