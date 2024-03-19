################# THIS IS THE EXAMPLE HOW TO USE RF MODEL IN NEW DATA #################
################# THIS IS THE EXAMPLE HOW TO USE RF MODEL IN NEW DATA #################
################# THIS IS THE EXAMPLE HOW TO USE RF MODEL IN NEW DATA #################
################# THIS IS THE EXAMPLE HOW TO USE RF MODEL IN NEW DATA #################
################# THIS IS THE EXAMPLE HOW TO USE RF MODEL IN NEW DATA #################
################# Prediction model and example files can be downloaded from https://zenodo.org/records/10258748 #################
library('e1071')
library('multiclassPairs')
library('openxlsx')
library('ggplot2')
library('survminer')
library('survival')

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
  
KMsurvival_plot<-ggsurvplot(fit,data=result_GSE181559_normal,pval = TRUE, #show p-value of log-rank test，
                              #conf.int = TRUE, #
                              pval.size=10,
                              legend.labs =  c( "Subgroup1",'Subgroup2','Subgroup3'),
                              palette = c( "#6688AB", "#97C17E",'#BA86B5'),
                              legend.title='', 
                              xlab = "Time in Day",   ###  customize X axis label.
                              #xlim=c(0,50),
                              break.x.by=1000, ###
                              ylab=paste0('Overall Survival'),
                              #ylab=paste0('Overall peak day'),
                              #ylab=paste0('SRAS−CoV−2 RNA +'),
                              surv.median.line = "hv", #
                              #palette = c( "blue","red"), ###  
                              #font.main = c(16, "bold", "darkblue"),
                              font.x = 40, # X 轴
                              font.y = 35, # c(14, "bold.italic", "darkred"),
                              font.tickslab = 30,# c(12, "plain", "darkgreen"), 
                              #conf.int.style = "step",  ###  customize style of confidence intervals,
                              risk.table = "abs_pct",  ###  absolute number and percentage at risk，
                              risk.table.y.text.col = T,###  colour risk table text annotations.
                              risk.table.y.text = F,###  show bars instead of names in text annotations in legend of risk table.
                              tables.y.text=F,
                              tables.x.text =F,
                              #risk.table.title="My title", ##  title name
                              fontsize=5, ## font size
                              ncensor.plot = F, #
                              #tables.theme=theme_cleantable(), # table
                              ggtheme = theme_classic()#
)
KMsurvival_plot$plot<-KMsurvival_plot$plot+
    theme(legend.text = element_text(size = 30),plot.margin = unit(c(2,2,0,2), "cm"),
          axis.text=element_text(colour="black"),
          axis.title=element_text(colour="black"))
  
KMsurvival_plot$table<-KMsurvival_plot$table+labs(x = NULL, y = NULL)+theme(axis.text=element_text(size=30,colour="black"),
                                                                              axis.title=element_text(size=30,colour="black"),
                                                                              
                                                                              legend.text = element_text(size = 25,colour="black"),
                                                                              plot.title = element_text(size=25),
                                                                              plot.margin = unit(c(0,2,2,2), "cm")) #  )
KMsurvival_plot

################# EXAMPLE END #################
################# EXAMPLE END #################
################# EXAMPLE END #################
################# EXAMPLE END #################
################# EXAMPLE END #################
