#### Figure 7 ####
library('e1071')
library(gplots)
library(ROCR)
library(multicore)
library(Rgtsp)
library(openxlsx)
library(multiclassPairs)
getwd()
load('/data/ylzhou/AMIS/Merged_Micor_RNA_NonNormalised.Rdata')
# load Merged_NonNormalised_matrix_RNA & MYCN_conditions_new_NonNormalised

MYCN_conditions_new_NonNormalised_train<-MYCN_conditions_new_NonNormalised[MYCN_conditions_new_NonNormalised$Splitgroup=='Train',]
MYCN_conditions_new_NonNormalised_test<-MYCN_conditions_new_NonNormalised[MYCN_conditions_new_NonNormalised$Splitgroup=='Test',]

Merged_NonNormalised_matrix_RNA_train<-Merged_NonNormalised_matrix_RNA[,MYCN_conditions_new_NonNormalised_train[,8]]
Merged_NonNormalised_matrix_RNA_test<-Merged_NonNormalised_matrix_RNA[,MYCN_conditions_new_NonNormalised_test[,8]]

object <- ReadData(Data = as.matrix(Merged_NonNormalised_matrix_RNA_train),
                   Labels =  MYCN_conditions_new_NonNormalised_train$Subgroup,
                   Platform = MYCN_conditions_new_NonNormalised_train$Platform_withRNA,
                   verbose = FALSE)
class(object) # "multiclassPairs_object"

filtered_genes <- filter_genes_TSP(data_object = object,
                                   filter = "one_vs_one",
                                   platform_wise = FALSE,
                                   featureNo = 1000,
                                   UpDown = TRUE,
                                   verbose = TRUE)
filtered_genes

# Let's train our model
classifier <- train_one_vs_rest_TSP(data_object = object,
                                    filtered_genes = filtered_genes,
                                    k_range = 5:50,
                                    include_pivot = FALSE,
                                    one_vs_one_scores = TRUE,
                                    platform_wise_scores = FALSE,
                                    seed = 1234,
                                    verbose = FALSE)
classifier

## Prediction
results_train <- predict_one_vs_rest_TSP(classifier = classifier,
                                         Data = object,
                                         tolerate_missed_genes = TRUE,
                                         weighted_votes = TRUE,
                                         classes = c("MYCN-Normal-Subgroup1",
                                         "MYCN-Normal-Subgroup2",
                                         "MYCN-Normal-Subgroup3"),
                                         verbose = TRUE)
# apply on the testing data
results_test <- predict_one_vs_rest_TSP(classifier = classifier,
                                        Data = as.matrix(Merged_NonNormalised_matrix_RNA_test),
                                        tolerate_missed_genes = TRUE,
                                        weighted_votes = TRUE,
                                        classes=c("MYCN-Normal-Subgroup1",
                                         "MYCN-Normal-Subgroup2",
                                         "MYCN-Normal-Subgroup3"),
                                        verbose = TRUE)

knitr::kable(head(results_test))

# Confusion Matrix and Statistics on training data
caret::confusionMatrix(data = factor(results_train$max_score, 
                                     levels = unique(object$data$Labels)),
                       reference = factor(object$data$Labels, 
                                          levels = unique(object$data$Labels)),
                       mode="everything")

results_test$max_score[1:10]
MYCN_conditions_new_NonNormalised_test$Subgroup[1:10]
colnames(Merged_NonNormalised_matrix_RNA_test)[1:10]
MYCN_conditions_new_NonNormalised_test$RNA[1:10]

# Confusion Matrix and Statistics on testing data
caret::confusionMatrix(data = factor(results_test$max_score, 
                                     levels = unique(object$data$Labels)),
                       reference = factor(MYCN_conditions_new_NonNormalised_test$Subgroup, 
                                          levels = unique(object$data$Labels)),
                       mode="everything")

#plot_binary_TSP(Data = object, # we are using the data object here
 #               classifier = classifier, 
  #              prediction = results_train, 
   #             classes = c("MYCN-Normal-Subgroup1",
    #                                     "MYCN-Normal-Subgroup2",
     #                                    "MYCN-Normal-Subgroup3"),
      #          #margin = c(0,5,0,10),
       #         title = "Training data")
1+1
# Random Forest scheme
# (500 trees here just for fast example)
genes_RF <- sort_genes_RF(data_object = object,
                          rank_data = TRUE,
                          platform_wise = TRUE,
                          num.trees = 10000, # more features, more tress are recommended
                          seed=123456, # for reproducibility
                          verbose = TRUE)

genes_RF # sorted genes object


# to get an idea of how many genes we will use
# and how many rules will be generated
summary_genes <- summary_genes_RF(sorted_genes_RF = genes_RF,
                                  genes_altogether = c(10,20,50,100,150,200,250,500,1000),
                                  genes_one_vs_rest = c(10,20,50,100,150,200,250,500,1000))
knitr::kable(summary_genes)

rules_RF <- sort_rules_RF(data_object = object, 
                          sorted_genes_RF = genes_RF,
                          genes_altogether = 200,
                          genes_one_vs_rest = 200, 
                          num.trees = 10000,# more rules, more tress are recommended 
                          seed=123456,
                          verbose = TRUE)
rules_RF # sorted rules object

# prepare the simple data.frame for the parameters I want to test
# names of arguments as column names
# this df has three sets (3 rows) of parameters
parameters <- data.frame(
  gene_repetition=c(200,100,50,20,10,5,1),
  rules_one_vs_rest=c(1000,1000,1000,1000,1000,1000,1000),
  rules_altogether=c(1000,1000,1000,1000,1000,1000,1000),
  run_boruta=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,TRUE), # I want to produce error in the 2nd trial
  plot_boruta = FALSE,
  num.trees=c(10000,10000,10000,10000,10000,10000,10000),
  stringsAsFactors = FALSE)

# parameters
# for overall and byclass possible options, check the help files
para_opt <- optimize_RF(data_object = object,
                        sorted_rules_RF = rules_RF,
                        parameters = parameters,
                        test_object = NULL,
                        overall = c("Accuracy","Kappa"), # wanted overall measurements 
                        byclass = c("F1"), # wanted measurements per class
                        verbose = TRUE)

para_opt # results object
# para_opt$summary # the df of with summarized information
knitr::kable(para_opt$summary)


# train the final model
# it is preferred to increase the number of trees and rules in case you have
# large number of samples and features
# for quick example, we have small number of trees and rules here
# based on the optimize_RF results we will select the parameters
RF_classifier <- train_RF(data_object = object,
                          sorted_rules_RF = rules_RF,
                          gene_repetition = 10,
                          rules_altogether = 1000,
                          rules_one_vs_rest = 1000,
                          run_boruta = TRUE, 
                          plot_boruta = FALSE,
                          probability = TRUE,
                          num.trees = 100000,
                          boruta_args = list(),
                          seed=123456,
                          verbose = TRUE)

#proximity_matrix_RF(object = object,
 #            classifier = RF_classifier, 
  #           plot = TRUE,
   #          return_matrix = FALSE, # if we need to get the matrix itself
    #         title = "Train",
     #        cluster_cols = TRUE)

# training accuracy
# get the prediction labels from the trained model
# if the classifier trained using probability   = FALSE
training_pred <- RF_classifier$RF_scheme$RF_classifier$predictions
if (is.factor(training_pred)) {
  x <- as.character(training_pred)
}

# if the classifier trained using probability   = TRUE
if (is.matrix(training_pred)) {
  x <- colnames(training_pred)[max.col(training_pred)]
}

# training accuracy
caret::confusionMatrix(data =factor(x),
                       reference = factor(object$data$Labels),
                       mode = "everything")

# apply on test data
results <- predict_RF(classifier = RF_classifier, 
                      Data = as.matrix(Merged_NonNormalised_matrix_RNA_test),
                      impute = TRUE) # can handle missed genes by imputation

# get the prediction labels
# if the classifier trained using probability   = FALSE
test_pred <- results$predictions
if (is.factor(test_pred)) {
  x <- as.character(test_pred)
}

# if the classifier trained using probability   = TRUE
if (is.matrix(test_pred)) {
  x <- colnames(test_pred)[max.col(test_pred)]
}

# training accuracy
caret::confusionMatrix(data = factor(x),
                       reference = factor(MYCN_conditions_new_NonNormalised_test$Subgroup),
                       mode = "everything") 

save( RF_classifier,file='RF_classifier.Rdata')

################# THIS IS THE EXAMPLE HOW TO USE RF MODEL IN NEW DATA #################
################# THIS IS THE EXAMPLE HOW TO USE RF MODEL IN NEW DATA #################
################# THIS IS THE EXAMPLE HOW TO USE RF MODEL IN NEW DATA #################
################# THIS IS THE EXAMPLE HOW TO USE RF MODEL IN NEW DATA #################
################# THIS IS THE EXAMPLE HOW TO USE RF MODEL IN NEW DATA #################
load('RF_classifier.Rdata')
GSE181559<-read.csv('/Users/yz3n18/Neuroblastoma_2022/Final_Figure0718/GSE181559_Neuroblastoma_tumor_totalRNA_TMM_FPKM.csv',
                      header = T,sep = ',',row.names = 1)
GSE181559<-GSE181559[,-1]
phe_GSE181559<-read.xlsx('/Users/yz3n18/Neuroblastoma_2022/Final_Figure0718/GSE181559.xlsx')
colnames(GSE181559)<-phe_GSE181559$X1
result_GSE181559<-predict_RF(classifier = RF_classifier,
                               Data=as.matrix(GSE181559),
                               impute = TRUE)
result_GSE181559<-(result_GSE181559$predictions)
if (is.matrix(result_GSE181559)){
    x<- colnames(result_GSE181559)[max.col(result_GSE181559)]
}
result_GSE181559<-as.data.frame(result_GSE181559)
result_GSE181559$Subgroup<-x
result_GSE181559$X1<-rownames(result_GSE181559)
result_GSE181559<-merge(result_GSE181559,phe_GSE181559,by='X1')
result_GSE181559_normal<-result_GSE181559[result_GSE181559$MYCN==0,]
fit<-survfit(Surv(Day,Event) ~ Subgroup, data=result_GSE181559_normal)
  
KMsurvival_plot<-ggsurvplot(fit,data=result_GSE181559_normal,pval = TRUE, #show p-value of log-rank test，显示log-rank分析得到的P值
                              #conf.int = TRUE, #添加置信
                              pval.size=10,
                              legend.labs =  c( "Subgroup1",'Subgroup2','Subgroup3'),
                              palette = c( "#6688AB", "#97C17E",'#BA86B5'),
                              legend.title='', 
                              xlab = "Time in Day",   ###  customize X axis label.自定义x的标time in years
                              #xlim=c(0,50),
                              break.x.by=1000, ###改变坐标轴间
                              ylab=paste0('Overall Survival'),
                              #ylab=paste0('Overall peak day'),
                              #ylab=paste0('SRAS−CoV−2 RNA +'),
                              surv.median.line = "hv", #添加中
                              #palette = c( "blue","red"), ###  自定义颜色
                              #font.main = c(16, "bold", "darkblue"),
                              font.x = 40, # X 轴
                              font.y = 35, # c(14, "bold.italic", "darkred"), y 轴特征
                              font.tickslab = 30,# c(12, "plain", "darkgreen"), 坐标轴大小
                              #conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样
                              risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)risk table
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
    theme(legend.text = element_text(size = 30),plot.margin = unit(c(2,2,0,2), "cm"),
          axis.text=element_text(colour="black"),
          axis.title=element_text(colour="black"))
  
KMsurvival_plot$table<-KMsurvival_plot$table+labs(x = NULL, y = NULL)+theme(axis.text=element_text(size=30,colour="black"),
                                                                              axis.title=element_text(size=30,colour="black"),
                                                                              
                                                                              legend.text = element_text(size = 25,colour="black"),
                                                                              plot.title = element_text(size=25),
                                                                              plot.margin = unit(c(0,2,2,2), "cm")) # 改变numer at risk 大小)
KMsurvival_plot

################# EXAMPLE END #################
################# EXAMPLE END #################
################# EXAMPLE END #################
################# EXAMPLE END #################
################# EXAMPLE END #################


# LOAD WESTERMAN ALK DATASET
Wester_ALK<-read.table('/data/ylzhou/AMIS/Westermann_ALK_logTPM.txt',header = T,sep='\t')
phe_Wester_ALK<-read.xlsx('/data/ylzhou/AMIS/8_Neuroblastoma_ALT_Westermann_144.xlsx')

colnames(Wester_ALK)
hist(Wester_ALK$fw2010nbd12)
1+1
Wester_ALK<-aggregate(x=Wester_ALK[,3:146],by=list(Wester_ALK$H.hugo),FUN=max)
rownames(Wester_ALK)<-Wester_ALK$Group.1
Wester_ALK<-Wester_ALK[,-1]
Wester_ALK_normal<-Wester_ALK[,tolower(phe_Wester_ALK[phe_Wester_ALK$MYCN=='normal',6])]


results_West <- predict_RF(classifier = RF_classifier, 
                      Data = as.matrix(Wester_ALK_normal),
                      impute = TRUE) # can handle missed genes by imputation
results_West
# if the classifier trained using probability   = FALSE
West_pred <- (results_West$predictions)
if (is.factor(West_pred)) {
  x <- as.character(West_pred)
}

# if the classifier trained using probability   = TRUE
if (is.matrix(West_pred)) {
  x <- colnames(West_pred)[max.col(West_pred)]
}
West_pred <- as.data.frame(results_West$predictions)
West_pred$Sample_2<-toupper(rownames(West_pred))
West_pred$Subgroup<-x
West_pred<-merge(West_pred,phe_Wester_ALK,by='Sample_2')
library(survival)
library(survminer)
fit<-survfit(Surv(Day,Event) ~Subgroup, data=West_pred)

KMsurvival_plot<-ggsurvplot(fit,data=West_pred,pval = TRUE, #show p-value of log-rank test，显示log-rank分析得到的P值
                            #conf.int = TRUE, #添加置信区间
                            pval.size=10,
                            legend.labs =  c( "Subgroup1",'Subgroup2','Subgroup3'),
                            palette = c( "#6688AB", "#97C17E",'#BA86B5'),
                            legend.title='', 
                            xlab = "Time in Day",   ###  customize X axis label.自定义x的标time in years
                            #xlim=c(0,50),
                            #break.x.by=1000, ###改变坐标轴间
                            ylab=paste0('Overall Survival'),
                            #ylab=paste0('Overall peak day'),
                            #ylab=paste0('SRAS−CoV−2 RNA +'),
                            surv.median.line = "hv", #添加中
                            #palette = c( "blue","red"), ###  自定义颜色
                            #font.main = c(16, "bold", "darkblue"),
                            font.x = 40, # X 轴
                            font.y = 35, # c(14, "bold.italic", "darkred"), y 轴特征
                            font.tickslab = 30,# c(12, "plain", "darkgreen"), 坐标轴大小
                            #conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)risk table
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
KMsurvival_plot
ggsave(KMsurvival_plot,file='Westerman_NonNormalized_KM_20230714.pdf',width=15,height = 12)
save(West_pred,file='West_pred_20230718.Rdata')

# TRAGET micro

TARGET_micro<-read.table('/data/ylzhou/AMIS/TARGET/Target_data_merging.txt',header=T,sep='\t',row.names = 1)
TARGET_micro_clinical<-read.table('/data/ylzhou/AMIS/TARGET/Target_data_merging_clinical.txt',header=T,sep='\t')


results_target_micro <- predict_RF(classifier = RF_classifier, 
                      Data = as.matrix(TARGET_micro),
                      impute = TRUE) # can handle missed genes by imputation
results_target_micro
# if the classifier trained using probability   = FALSE
Target_micro_pred <- (results_target_micro$predictions)


# if the classifier trained using probability   = TRUE
if (is.matrix(Target_micro_pred)) {
  x <- colnames(Target_micro_pred)[max.col(Target_micro_pred)]
}

Target_micro_pred<-as.data.frame(Target_micro_pred)

Target_micro_pred$Subgroup<-x
Target_micro_pred$TARGET.USI<-rownames(Target_micro_pred)

Target_micro_pred<-merge(Target_micro_pred,TARGET_micro_clinical,by='TARGET.USI')
colnames(Target_micro_pred)

fit<-survfit(Surv(Overall.Survival.Time.in.Days,Vital.Status.1) ~Subgroup, data=Target_micro_pred)

KMsurvival_plot<-ggsurvplot(fit,data=Target_micro_pred,pval = TRUE, #show p-value of log-rank test，显示log-rank分析得到的P值
                            #conf.int = TRUE, #添加置信区间
                            pval.size=10,
                            legend.labs =  c( "Subgroup1",'Subgroup2','Subgroup3'),
                            palette = c( "#6688AB", "#97C17E",'#BA86B5'),
                            legend.title='', 
                            xlab = "Time in Day",   ###  customize X axis label.自定义x的标time in years
                            #xlim=c(0,50),
                            #break.x.by=1000, ###改变坐标轴间
                            ylab=paste0('Overall Survival'),
                            #ylab=paste0('Overall peak day'),
                            #ylab=paste0('SRAS−CoV−2 RNA +'),
                            surv.median.line = "hv", #添加中
                            #palette = c( "blue","red"), ###  自定义颜色
                            #font.main = c(16, "bold", "darkblue"),
                            font.x = 40, # X 轴
                            font.y = 35, # c(14, "bold.italic", "darkred"), y 轴特征
                            font.tickslab = 30,# c(12, "plain", "darkgreen"), 坐标轴大小
                            #conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)risk table
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
KMsurvival_plot

save(Target_micro_pred,file='Target_micro_pred_20230718.Rdata')
# TRAGET RNA

TARGET_RNA<-read.table('/data/ylzhou/AMIS/TARGET/Target_RNA-seq_merging_2.txt',header=T,sep='\t')
TARGET_rna_clinical<-read.table('/data/ylzhou/AMIS/TARGET/Target_Clinical.csv',header=T,sep=',')
TARGET_sampleID<-read.xlsx('/data/ylzhou/AMIS/TARGET/RNA_Sample_ID.xlsx')
TARGET_rna_clinical<-merge(TARGET_rna_clinical,TARGET_sampleID,by='Sample_2')


colnames(TARGET_RNA)
TARGET_RNA<-aggregate(x=TARGET_RNA[,2:156],by=list(TARGET_RNA$Gene),FUN=max)
rownames(TARGET_RNA)<-TARGET_RNA$Group.1
TARGET_RNA<-TARGET_RNA[,-1]


results_target_rna <- predict_RF(classifier = RF_classifier, 
                      Data = as.matrix(TARGET_RNA),
                      impute = TRUE) # can handle missed genes by imputation

# if the classifier trained using probability   = FALSE
Target_rna_pred <- (results_target_rna$predictions)


# if the classifier trained using probability   = TRUE
if (is.matrix(Target_rna_pred)) {
  x <- colnames(Target_rna_pred)[max.col(Target_rna_pred)]
}

Target_rna_pred<-as.data.frame(Target_rna_pred)

Target_rna_pred$Subgroup<-x
Target_rna_pred$Sample<-rownames(Target_rna_pred)
Target_rna_pred<-merge(Target_rna_pred,TARGET_rna_clinical,by='Sample')
colnames(Target_rna_pred)


fit<-survfit(Surv(Overall.Survival.Time.in.Days,Vital.Status) ~Subgroup, data=Target_rna_pred)

KMsurvival_plot<-ggsurvplot(fit,data=Target_rna_pred,pval = TRUE, #show p-value of log-rank test，显示log-rank分析得到的P值
                            #conf.int = TRUE, #添加置信区间
                            pval.size=10,
                            legend.labs =  c( "Subgroup1",'Subgroup2','Subgroup3'),
                            palette = c( "#6688AB", "#97C17E",'#BA86B5'),
                            legend.title='', 
                            xlab = "Time in Day",   ###  customize X axis label.自定义x的标time in years
                            #xlim=c(0,50),
                            #break.x.by=1000, ###改变坐标轴间
                            ylab=paste0('Overall Survival'),
                            #ylab=paste0('Overall peak day'),
                            #ylab=paste0('SRAS−CoV−2 RNA +'),
                            surv.median.line = "hv", #添加中
                            #palette = c( "blue","red"), ###  自定义颜色
                            #font.main = c(16, "bold", "darkblue"),
                            font.x = 40, # X 轴
                            font.y = 35, # c(14, "bold.italic", "darkred"), y 轴特征
                            font.tickslab = 30,# c(12, "plain", "darkgreen"), 坐标轴大小
                            #conf.int.style = "step",  ###  customize style of confidence intervals,改变置信区间的样
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk，这里以n(%)risk table
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
KMsurvival_plot
save(Target_rna_pred,file='Target_rna_pred_20230718.Rdata')


# GSE49710
load('/data/ylzhou/AMIS/GSE49710_NonNormalised.Rdata')
phe_GSE49710<-read.csv('/data/ylzhou/AMIS/phe_GSE49710.csv',header=T,sep=',')
results_gse49710 <- predict_RF(classifier = RF_classifier, 
                      Data = as.matrix(GSE49710),
                      impute = TRUE) # can handle missed genes by imputation

# if the classifier trained using probability   = FALSE
GSE49710_pred <- (results_gse49710$predictions)


# if the classifier trained using probability   = TRUE
if (is.matrix(GSE49710_pred)) {
  x <- colnames(GSE49710_pred)[max.col(GSE49710_pred)]
}

GSE49710_pred<-as.data.frame(GSE49710_pred)

rownames(GSE49710_pred)<-gsub('gProcessedSignal.','',rownames(GSE49710_pred))

rownames(GSE49710_pred)<-sapply(rownames(GSE49710_pred),function(x){
  strsplit(x,'_')[[1]][3]
})
GSE49710_pred$RNA<-rownames(GSE49710_pred)
GSE49710_pred$Subgroup<-x

GSE49710_pred<-merge(GSE49710_pred,MYCN_conditions_new_NonNormalised,by='RNA')
colnames(GSE49710_pred)
table(GSE49710_pred$Subgroup.x,GSE49710_pred$Subgroup.y)

table(GSE49710_pred$MYCN_conditions,GSE49710_pred$Subgroup.y)
save(GSE49710_pred,file='GSE49710_pred_0718.Rdata')

