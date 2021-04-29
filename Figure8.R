library(caret)
BiocManager::install('glmnet')
library(glmnet)
library(superml)


#### 43 Predictors ####

dataObj1 = datasetObject_sil_test$sillhouette$GSE120572 # 394 samples
dataObj2 = datasetObject_sil_test$sillhouette$GSE112447 # 57 samples
dataObj3 = datasetObject_sil_test$sillhouette$GSE73517 # 105 samples 
dataObj5 = datasetObject_sil_test$sillhouette$GSE49710 # 498 samples
dataObj7 = datasetObject_sil_test$sillhouette$GSE45547 # 649 samples
dataObj8 = datasetObject_sil_test$sillhouette$GSE19274 #100 samples
dataObj9 = datasetObject_sil_test$sillhouette$GSE16237 # 51 samples
dataObj10 = datasetObject_sil_test$sillhouette$GSE14880 # 29 samples Delete it for subgroup2 as less than 5 sample are subgroup2 in this data
dataObj11 = datasetObject_sil_test$sillhouette$GSE13136 #30 samples
dataObj13 = datasetObject_sil_test$sillhouette$GSE27608 # 47 samples
dataObj14 = datasetObject_sil_test$sillhouette$GSE21713 # 40 samples
dataObj15 = datasetObject_sil_test$sillhouette$GSE85047 #283 samples
dataObj16 = datasetObject_sil_test$sillhouette$GSE16476 #88 samples
dataObj17 = datasetObject_sil_test$sillhouette$GSE79910 #27 samples Delete it for subgroup2
dataObj18 = datasetObject_sil_test$sillhouette$E_MTAB_1781 #683 samples
dataObj19 = datasetObject_sil_test$sillhouette$E_MTAB_8248 #222 samples
dataObj20 = datasetObject_sil_test$sillhouette$E_MTAB_179 #477 samples
dataObj21 = datasetObject_sil_test$sillhouette$E_MEXP_3517 #42 samples

validation_datasets <- list(dataObj1, dataObj2,dataObj3,
                            dataObj5,dataObj7,dataObj8,
                            dataObj9,dataObj11,dataObj10,
                            dataObj13,dataObj14,dataObj15,
                            dataObj16,dataObj18,dataObj17,
                            dataObj19,dataObj20,dataObj21)
names(validation_datasets)<-c(dataObj1$formattedName, dataObj2$formattedName,dataObj3$formattedName,
                              dataObj5$formattedName,dataObj7$formattedName,dataObj8$formattedName,
                              dataObj9$formattedName,dataObj11$formattedName,dataObj10$formattedName,
                              dataObj13$formattedName,dataObj14$formattedName,dataObj15$formattedName,
                              dataObj16$formattedName,dataObj18$formattedName,dataObj17$formattedName,
                              dataObj19$formattedName,dataObj20$formattedName,dataObj21$formattedName)


exampleMetaObj=list() 
exampleMetaObj$originalData <- discovery_datasets
exampleMetaObj_validation<-list()
exampleMetaObj_validation$originalData <- validation_datasets

checkDataObject(exampleMetaObj, "Meta", "Pre-Analysis") ###Check your metaObject before MetaAnalysis using 
checkDataObject(exampleMetaObj_validation, "Meta", "Pre-Analysis")
exampleMetaObj <- runMetaAnalysis(exampleMetaObj, maxCores=4,runLeaveOneOutAnalysis = FALSE)
exampleMetaObj <- filterGenes(exampleMetaObj, isLeaveOneOut = F, FDRThresh = 0.05,effectSizeThresh = 1.3)
exampleMetaObj$filterResults$FDR0.05_es1.3_nStudies1_looaFALSE_hetero0$posGeneNames<-c(intersect(rownames(test_1), 
                                                                                                 exampleMetaObj$filterResults$FDR0.05_es1.3_nStudies1_looaFALSE_hetero0$posGeneNames))
exampleMetaObj$filterResults$FDR0.05_es1.3_nStudies1_looaFALSE_hetero0$negGeneNames<-c(intersect(rownames(test_1), 
                                                                                                 exampleMetaObj$filterResults$FDR0.05_es1.3_nStudies1_looaFALSE_hetero0$negGeneNames))

forwardRes <- forwardSearch( metaObject = exampleMetaObj, 
                             filterObject = exampleMetaObj$filterResults$FDR0.05_es1.3_nStudies1_looaFALSE_hetero0,
                             forwardThresh = 0.0000000000001)

#### Supplementary Figure 6a ####
summaryROCPlot(metaObject = exampleMetaObj_validation, 
               filterObject = exampleMetaObj$filterResults$FDR0.05_es1.3_nStudies1_looaFALSE_hetero0, 
               bootstrapReps = 100)

#### SVM ####

merge_data_batch_MYCN_amp_train_train1<-data.frame(train_1,train_2,train_3)
merge_data_batch_MYCN_amp_test_test1<-data.frame(test_1,test_2,test_3)
merge_data_batch_MYCN_amp_train_train1_sub1pos<-t(merge_data_batch_MYCN_amp_train_train1[c("SCN2A","LONRF2","FRS3","CPEB4","SNAP25","PMP22","OLA1","C14orf132"),]) ### 
merge_data_batch_MYCN_amp_train_train1_sub1neg<-t(merge_data_batch_MYCN_amp_train_train1[c("TMEM109","LDHA","NDE1","PIM2","MRPL11","TNFRSF10B","COLEC12","TAF10","ELL","SIVA1"),])
merge_data_batch_MYCN_amp_train_train1_sub2pos<-t(merge_data_batch_MYCN_amp_train_train1[c("KIF4A","COX8A","UHRF1","ODC1","SNRPD1","CDCA4","SKP2","HMGB3","FANCI","CDK2"),])
merge_data_batch_MYCN_amp_train_train1_sub2neg<-t(merge_data_batch_MYCN_amp_train_train1[c("PLA2G4C","NDEL1","ANXA2","PLSCR3"),])
merge_data_batch_MYCN_amp_train_train1_sub3pos<-t(merge_data_batch_MYCN_amp_train_train1[c("IFITM2","ANXA11","STAT5A","FYCO1","LMNA","RARRES3","PLSCR4","SLC10A3","TMCO4","AGA"),])
merge_data_batch_MYCN_amp_train_train1_sub3neg<-t(merge_data_batch_MYCN_amp_train_train1[c("KBTBD7"),])

merge_data_batch_MYCN_amp_test_test1_sub1pos<-t(merge_data_batch_MYCN_amp_test_test1[c("SCN2A","LONRF2","FRS3","CPEB4","SNAP25","PMP22","OLA1","C14orf132"),])
merge_data_batch_MYCN_amp_test_test1_sub1neg<-t(merge_data_batch_MYCN_amp_test_test1[c("TMEM109","LDHA","NDE1","PIM2","MRPL11","TNFRSF10B","COLEC12","TAF10","ELL","SIVA1"),])
merge_data_batch_MYCN_amp_test_test1_sub2pos<-t(merge_data_batch_MYCN_amp_test_test1[c("KIF4A","COX8A","UHRF1","ODC1","SNRPD1","CDCA4","SKP2","HMGB3","FANCI","CDK2"),])
merge_data_batch_MYCN_amp_test_test1_sub2neg<-t(merge_data_batch_MYCN_amp_test_test1[c("PLA2G4C","NDEL1","ANXA2","PLSCR3"),])
merge_data_batch_MYCN_amp_test_test1_sub3pos<-t(merge_data_batch_MYCN_amp_test_test1[c("IFITM2","ANXA11","STAT5A","FYCO1","LMNA","RARRES3","PLSCR4","SLC10A3","TMCO4","AGA"),])
merge_data_batch_MYCN_amp_test_test1_sub3neg<-t(merge_data_batch_MYCN_amp_test_test1[c("KBTBD7"),])
#### Geometric Mean ####
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}
train_GM<-matrix(nrow=length(merge_data_batch_MYCN_amp_train_train1[1,]),# nrow= number of sample,ncol=6 GM values+Group
                 ncol=7,data=NA)
colnames(train_GM)<-c('Sub1Pos','Sub1Neg',
                      'Sub2Pos','Sub2Neg',
                      'Sub3Pos','Sub3Neg','Group')

dim(train_3)
train_GM[,1]<-apply(merge_data_batch_MYCN_amp_train_train1_sub1pos,1,gm_mean) ## add GM of subgoup1 postive
train_GM[,2]<-apply(merge_data_batch_MYCN_amp_train_train1_sub1neg,1,gm_mean) ## add GM of subgoup1 negative
train_GM[,3]<-apply(merge_data_batch_MYCN_amp_train_train1_sub2pos,1,gm_mean)
train_GM[,4]<-apply(merge_data_batch_MYCN_amp_train_train1_sub2neg,1,gm_mean)
train_GM[,5]<-apply(merge_data_batch_MYCN_amp_train_train1_sub3pos,1,gm_mean)
train_GM[,6]<-apply(merge_data_batch_MYCN_amp_train_train1_sub3neg,1,gm_mean)
train_GM[,7]<-c(rep(1,1034),rep(2,316),rep(3,370))### add group information silhoutte 
train_GM[,7]<-c(rep(1,1034),rep(2,559),rep(3,567))### add group information All
rownames(train_GM)<-rownames(merge_data_batch_MYCN_amp_train_train1_sub1neg)


test_GM<-matrix(nrow=length(merge_data_batch_MYCN_amp_test_test1[1,]),
                ncol=7,data=NA)
colnames(test_GM)<-c('Sub1Pos','Sub1Neg',
                     'Sub2Pos','Sub2Neg',
                     'Sub3Pos','Sub3Neg','Group')
dim(test_3)
test_GM[,1]<-apply(merge_data_batch_MYCN_amp_test_test1_sub1pos,1,gm_mean) ## add GM of subgoup1 postive
test_GM[,2]<-apply(merge_data_batch_MYCN_amp_test_test1_sub1neg,1,gm_mean) ## add GM of subgoup1 negative
test_GM[,3]<-apply(merge_data_batch_MYCN_amp_test_test1_sub2pos,1,gm_mean)
test_GM[,4]<-apply(merge_data_batch_MYCN_amp_test_test1_sub2neg,1,gm_mean)
test_GM[,5]<-apply(merge_data_batch_MYCN_amp_test_test1_sub3pos,1,gm_mean)
test_GM[,6]<-apply(merge_data_batch_MYCN_amp_test_test1_sub3neg,1,gm_mean)
test_GM[,7]<-c(rep(1,378),rep(2,108),rep(3,164))### add group information
test_GM[,7]<-c(rep(1,378),rep(2,268),rep(3,279))### add group information
rownames(test_GM)<-rownames(merge_data_batch_MYCN_amp_test_test1_sub1neg)
train_GM<-as.data.frame(train_GM)
test_GM<-as.data.frame(test_GM) ## 

#### calculate Z-score ####
train_GM_zscore<-train_GM
test_GM_zscore<-test_GM

zscore <- function(x){
  (x-mean(x))/sd(x)
}
train_GM_zscore[,1:6]<-apply(train_GM_zscore[,1:6], 2, zscore)# 
test_GM_zscore[,1:6]<-apply(test_GM_zscore[,1:6], 2, zscore)

#### SVM ####
train_zscore.x<- model.matrix(Group~., train_GM_zscore)[,-1]
train_zscore.y <- as.factor(train_GM_zscore$Group) ### 
# Convert the outcome (class) to a numerical variable
set.seed(123)
svm_model_radial_tune_zscore_removesomegene<-tune(svm, train.x = train_zscore.x,train.y=train_zscore.y, kernel = "radial",
                                                  ranges = list(gamma = seq(0.99,1.07,length.out=200), cost =seq(45,49,length.out=200)))
summary(svm_model_radial_tune_zscore_removesomegene)
# gamma 1.040526; cost 47.15789
svm_model_radial_after_tune_zscore_removesomegene <- svm(train_zscore.x,train_zscore.y, kernel="radial", cost=46.78894, 
                                                         gamma=1.043065) 
pred_zscore_removesomegene <- predict(svm_model_radial_after_tune_zscore_removesomegene,
                                      train_zscore.x) ###
table(pred_zscore_removesomegene,train_zscore.y)
test_zscore.x <- model.matrix(Group ~., test_GM_zscore)[,-1]
test_zscore.y<-as.factor(test_GM_zscore$Group)
pred_zscore_removesomegene<- predict(svm_model_radial_after_tune_zscore_removesomegene,test_zscore.x)
table(pred_zscore_removesomegene,test_zscore.y)
svm_model_radial_after_tune_zscore_removesomegene <- svm(train_zscore.x,train_zscore.y, kernel="radial", cost=46.78894, 
                                                         gamma=1.043065,probability=T) 
pred_zscore_removesomegene_score<- predict(svm_model_radial_after_tune_zscore_removesomegene_score,test_zscore.x,
                                           probability=T)


#### Figure 8a ####
library(ggplot2)
library(ggExtra)
Probability_Test<-as.data.frame(attr(pred_zscore_removesomegene_score, "probabilities")) # Test
Probability_Test$Sample<-rownames(test_GM)
Probability_Test$Group<-pred_zscore_removesomegene
Probability_Test$Cluster<-test_GM_zscore$Group
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
ggsave(file=paste0('Test_Probability2.pdf'),width = 6,height = 6,print(p1),
       onefile=FALSE)


#### XGBoost ####
# Random search
BiocManager::install('xgboost')
library(xgboost)

#preparing matrix 
dim(train.x)
train_XGB_labels <- as.numeric(train_GM_zscore$Group)-1
test_XGB_label <- as.numeric(test_GM_zscore$Group)-1
train.x<-model.matrix(Group~., train_GM_zscore)[,-1]
test.x<-model.matrix(Group~., test_GM_zscore)[,-1]

dtrain <- xgb.DMatrix(data = train.x,label = train_XGB_labels) 
dtest <- xgb.DMatrix(data = test.x,label=test_XGB_label)
watchlist <- list(train=dtrain, test=dtest)
#default parameters
params <- list(booster = "gbtree", objective = "multi:softmax", 
               eta=0.3, gamma=0, max_depth=10, min_child_weight=1, nthread = 4,
               subsample=1, colsample_bytree=1,num_class = 3)

# Take start time to measure time of random search algorithm
start.time <- Sys.time()
# Create empty lists
lowest_error_list = list()
parameters_list = list()

?runif
# Create 10000 rows with random hyperparameters
set.seed(123)
for (iter in 1:10000){
  param <- list(booster = "gbtree",
                objective = "multi:softmax",
                max_depth = sample(3:10, 1), #maximum tree depth
                eta = runif(1, .01, .4), #?runif  0.01 to 0.3  learning rate
                gamma=runif(1, .01, .4),# minimum loss reduction required to introduce a split in the classifier tree
                subsample = runif(1, .5, 1),
                colsample_bytree = runif(1, .5, 1),
                min_child_weight = sample(0:10, 1), #minimum child weight
                max_delta_step=runif(1, .6, 1),  #Maximum delta step
                lambda=runif(1, .01, .4),# L2 PENALTY
                nthread=6,
                tree_method=sample(c('exact','approx'),1),
                seed=123
  )
  parameters <- as.data.frame(param)
  parameters_list[[iter]] <- parameters
}

# Create object that contains all randomly created hyperparameters
parameters_df = do.call(rbind, parameters_list)

start.time <- Sys.time()
# Use randomly created parameters to create 10,000 XGBoost-models
for (row in 1:nrow(parameters_df)){
  set.seed(123)
  mdcv <- xgb.train(data=dtrain,
                    booster = "gbtree",
                    objective = "multi:softmax",
                    max_depth = parameters_df$max_depth[row],
                    eta = parameters_df$eta[row],
                    gamma=parameters_df$gamma[row],
                    subsample = parameters_df$subsample[row],
                    colsample_bytree = parameters_df$colsample_bytree[row],
                    min_child_weight = parameters_df$min_child_weight[row],
                    max_delta_step=parameters_df$max_delta_step[row],
                    lambda=parameters_df$lambda[row],
                    tree_method=parameters_df$tree_method[row],
                    nthread=6,
                    nrounds= 5000,
                    early_stopping_rounds= 50,
                    print_every_n = 100,
                    nfold = 10,
                    num_class = 3,
                    watchlist = list(train= dtrain, val= dtest)
  )
  lowest_error <- as.data.frame(1 - min(mdcv$evaluation_log$val_mlogloss))
  lowest_error_list[[row]] <- lowest_error
}
lowest_error_df = do.call(rbind, lowest_error_list)
# Bind columns of accuracy values and random hyperparameter values
randomsearch = cbind(lowest_error_df, parameters_df)
# Quickly display highest accuracy
max(randomsearch$`1 - min(mdcv$evaluation_log$val_mlogloss)`)

# Stop time and calculate difference
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
randomsearch_copy<-randomsearch
library(dplyr)
# Prepare table
randomsearch <- as.data.frame(randomsearch) %>%
  rename(val_acc = `1 - min(mdcv$evaluation_log$val_mlogloss)`) %>%
  arrange(-val_acc)
##
set.seed(123)
params <- list(booster = "gbtree", 
               objective = "multi:softmax", # multi:softprob,softmax
               max_depth = randomsearch[1,]$max_depth,
               eta = randomsearch[1,]$eta,
               subsample = randomsearch[1,]$subsample,
               colsample_bytree = randomsearch[1,]$colsample_bytree,
               min_child_weight = randomsearch[1,]$min_child_weight,
               gamma=randomsearch[1,]$gamma,
               max_delta_step=randomsearch[1,]$max_delta_step,
               lambda=randomsearch[1,]$lambda,
               tree_method=randomsearch[1,]$tree_method
)
xgb_tuned <- xgb.train(params = params,
                       data = dtrain,
                       nrounds =5000,
                       print_every_n = 100,
                       early_stopping_rounds = 50,
                       nfold = 10,
                       num_class = 3,nthread=6,watchlist = list(train= dtrain, val= dtest))

# Make prediction on dvalid
test_GM_XGBOOST<-test_GM_zscore
test_GM_XGBOOST$pred_survived_tuneProb <- predict(xgb_tuned, dtest)
test_GM_XGBOOST_Prob<-predict(xgb_tuned, dtest)

#### Comparison with public subgroups ####
# GSE85047 Califano Cancer Discov #
library(openxlsx)
library(sankeywheel)
getwd()

GSE85047_Califano_sankeywheel<-read.xlsx('/Comparison_Public_Subgroup/Andrea_Califano_Cancer_Discovery_TARGET_GSE85047.xlsx',
                                       sheet = 5)
GSE49710_Westermann_sankeywheel<-read.xlsx('Comparison_Public_Subgroup/Westermann_Nature_Cancer_GSE49710.xlsx',
                                           sheet = 4)
#### Figure 8c & Supplementary Figure 6b ####
sankeywheel(from = GSE85047_Califano_sankeywheel$fro,
            to = GSE85047_Califano_sankeywheel$to,
            weight =GSE85047_Califano_sankeywheel$weight,
            type = "sankey", 
            width = "100%",
            title = "")
#### Figure 8d ####

{GSE49710_Westermann_HR_sub3<-GSE49710_Westermann[GSE49710_Westermann$Signature=='MNA-HR',]
  GSE49710_Westermann_HR_sub3<-GSE49710_Westermann_mes_sub3[GSE49710_Westermann_mes_sub3$Zhou_Subgroup=='Subgroup3',]
  GSE49710_Westermann_HR_sub1<-GSE49710_Westermann[GSE49710_Westermann$Signature=='MNA-HR',]
  GSE49710_Westermann_HR_sub1<-GSE49710_Westermann_HR_sub1[GSE49710_Westermann_HR_sub1$Zhou_Subgroup=='Subgroup1',]
  GSE49710_Westermann_HR_sub2<-GSE49710_Westermann[GSE49710_Westermann$Signature=='MNA-HR',]
  GSE49710_Westermann_HR_sub2<-GSE49710_Westermann_HR_sub2[GSE49710_Westermann_HR_sub2$Zhou_Subgroup=='Subgroup2',]
  GSE49710_Westermann_HR<-rbind(GSE49710_Westermann_HR_sub3,
                                    GSE49710_Westermann_HR_sub1,
                                    GSE49710_Westermann_HR_sub2)
}
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
fit<-survfit(Surv(Day_survival,Event) ~ Zhou_Subgroup, data=GSE49710_Westermann_HR)
fit<-survfit(Surv(Day_survival,Event) ~ Zhou_Subgroup, data=GSE49710_Westermann_LR)

KMsurvival_plot<-ggsurvplot(fit,data=GSE49710_Westermann_HR,pval = TRUE, #show p-value of log-rank test
                            #conf.int = TRUE, #
                            pval.size=10,
                            legend.labs =  c("Subgroup1",'Subgroup2','Subgroup3'),
                            legend.title='', 
                            xlab = "Time in Day",   ###  customize X axis label
                            #xlim=c(0,50),
                            break.x.by=1000, ###
                            ylab=paste0('Overall survival'),
                            surv.median.line = "hv", 
                            palette = c('#7CAE00','#00BFC4', '#C77CFF'), 
                            #font.main = c(16, "bold", "darkblue"),
                            font.x = 40, # X 
                            font.y = 35, # c(14, "bold.italic", "darkred"), 
                            font.tickslab = 30,# c(12, "plain", "darkgreen"), 
                            #conf.int.style = "step",  ###  customize style of confidence intervals
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk
                            risk.table.y.text.col = T,###  colour risk table text annotations.
                            risk.table.y.text = F,###  show bars instead of names in text annotations in legend of risk table
                            tables.y.text=F,
                            tables.x.text =F,
                            #risk.table.title="My title", ## 
                            fontsize=5, ##
                            ncensor.plot = F, 
                            #tables.theme=theme_cleantable(), # 
                            ggtheme = theme_classic()#
)
KMsurvival_plot$plot<-KMsurvival_plot$plot+
  theme(legend.text = element_text(size = 30),plot.margin = unit(c(2,2,0,2), "cm"))

KMsurvival_plot$table<-KMsurvival_plot$table+labs(x = NULL, y = NULL)+theme(axis.text=element_text(size=30),
                                                                            axis.title=element_text(size=30),
                                                                            
                                                                            legend.text = element_text(size = 25),
                                                                            plot.title = element_text(size=25),
                                                                            plot.margin = unit(c(0,2,2,2), "cm"))

ggsave(file=paste0('Figure8d.pdf'),width = 15,height = 12,print(KMsurvival_plot),
       onefile=FALSE)
