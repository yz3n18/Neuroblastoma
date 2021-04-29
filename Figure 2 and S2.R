#### Figure 2 ####
#### Figure 2c/d & Supplementary Figure 2c  ####
BiocManager::install('survminer')
library(survminer)
library(survival)
survival_data<-read.csv('/Users/yihuawang/Neuroblastoma/Hu/Survival_ALL.csv')
rownames(survival_data)<-survival_data$Sample
survival_data<-na.omit(survival_data)
survival_data$Group<-NA
survival_data[colnames(train_1_sil),4]<-'Subgroup1' 
survival_data[colnames(train_2_sil),4]<-'Subgroup2' 
survival_data[colnames(train_3_sil),4]<-'Subgroup3' 
survival_data[colnames(merge_data_batch_MYCN_amp_train),4]<-'MYCN-AMP' 
survival_data[colnames(test_1_sil),4]<-'Subgroup1'
survival_data[colnames(test_2_sil),4]<-'Subgroup2' 
survival_data[colnames(test_3_sil),4]<-'Subgroup3' 
survival_data[colnames(merge_data_batch_MYCN_amp_test),4]<-'MYCN-AMP' 
survival_data<-na.omit(survival_data)
#colnames(survival_data)
#survival_data_AMP<-survival_data[survival_data[,4]!='MYCN-control',]
#survival_data_AMP<-survival_data[survival_data[,4]!='MYCN-AMP',]
#survival_data_AMP<-survival_data

fit<-survfit(Surv(Day_Survival,Event) ~ Group, data=survival_data)
KMsurvival_plot<-ggsurvplot(fit,data=survival_data,pval = TRUE, #show p-value of log-rank test
                            #conf.int = TRUE, #
                            pval.size=10,
                            legend.labs =  c("MYCN-AMP","Subgroup1",'Subgroup2','Subgroup3'),
                            legend.title='', 
                            xlab = "Time in Day",   ###  customize X axis label
                            #xlim=c(0,50),
                            break.x.by=1000, ###
                            ylab=paste0('Overall survival'),
                            #ylab=paste0('Overall peak day')
                            surv.median.line = "hv", #add line of median survival day
                            #palette = c("#7CAE00", "#00BFC4", "#C77CFF"), ###  
                            #font.main = c(16, "bold", "darkblue"),
                            font.x = 40, # X 
                            font.y = 35, # c(14, "bold.italic", "darkred"), y 
                            font.tickslab = 30,# c(12, "plain", "darkgreen"), 
                            #conf.int.style = "step",  ###  customize style of confidence intervals
                            risk.table = "abs_pct",  ###  absolute number and percentage at risk
                            risk.table.y.text.col = T,###  colour risk table text annotations.
                            risk.table.y.text = F,###  show bars instead of names in text annotations in legend of risk table.
                            tables.y.text=F,
                            tables.x.text =F,
                            #risk.table.title="My title", ## 
                            fontsize=5, ## table font
                            ncensor.plot = F, #
                            #tables.theme=theme_cleantable(), 
                            ggtheme = theme_classic()#
)
KMsurvival_plot$plot<-KMsurvival_plot$plot+
  theme(legend.text = element_text(size = 30),plot.margin = unit(c(2,2,0,2), "cm"))

KMsurvival_plot$table<-KMsurvival_plot$table+labs(x = NULL, y = NULL)+theme(axis.text=element_text(size=30),
                                                                            axis.title=element_text(size=30),
                                                                            
                                                                            legend.text = element_text(size = 25),
                                                                            plot.title = element_text(size=25),
                                                                            plot.margin = unit(c(0,2,2,2), "cm")))

ggsave(file=paste0('All_silhouette_Km_3_survival_Train+Test.pdf'),width = 15,height = 12,print(KMsurvival_plot),
       onefile=FALSE)



#### Figure 2e ####
survival_data_AMP_Cox_AURKA_MYCN_all<-read.csv('Survival/Cox/survival_data_AMP_Cox_AURKA_MYC_all.csv',
                                               header = T,sep = ',')

survival_data_AMP_Cox_AURKA_MYCN_Risk<-survival_data_AMP_Cox_AURKA_MYCN_all
survival_data_AMP_Cox_AURKA_MYCN_Risk$Gender <- factor(survival_data_AMP_Cox_AURKA_MYCN_Risk$Gender,levels = c('Male','Female'))
survival_data_AMP_Cox_AURKA_MYCN_Risk$High_risk <- factor(survival_data_AMP_Cox_AURKA_MYCN_Risk$High_risk,levels = c('Low risk','High risk'))
survival_data_AMP_Cox_AURKA_MYCN_Risk$Subgroup <- factor(survival_data_AMP_Cox_AURKA_MYCN_Risk$Subgroup,levels = c('Subgroup1','Subgroup2','Subgroup3'))
survival_data_AMP_Cox_AURKA_MYCN_Risk$INSS_Stage <- factor(survival_data_AMP_Cox_AURKA_MYCN_Risk$INSS_Stage,
                                                           levels = c('Stage2','Stage3','Stage4','Stage4S','Stage1'))

Basurv<-Surv(time=survival_data_AMP_Cox_AURKA_MYCN_Risk$Day_Survival,event = survival_data_AMP_Cox_AURKA_MYCN_Risk$Event)
Gcox_Risk<-coxph(Basurv~High_risk+Subgroup,data = survival_data_AMP_Cox_AURKA_MYCN_Risk)

pdf('Figure_2e.pdf',width=15,height = 6)
ggforest(Gcox_Risk,data = survival_data_AMP_Cox_AURKA_MYCN_Risk,cpositions = c(0.02,0.22,0.4),main = 'Hazard ratio',
         fontsize =1.8,refLabel = 'reference',noDigits = 2)
dev.off()

#### Figure 2f & Supplementary Figure 2d ####
library(pec)
library(prodlim)
library(survival)
Models <- list("Subgroup"=coxph(Basurv~Subgroup,x=TRUE,y=TRUE,
                                data = survival_data_AMP_Cox_AURKA_MYCN_all),
               "Risk"=coxph(Basurv~High_risk,x=TRUE,y=TRUE,
                            data = survival_data_AMP_Cox_AURKA_MYCN_all)
)

PredError <- pec(object=Models,
                 formula=Surv(Day_Survival,Event)~Subgroup+High_risk,
                 data=survival_data_AMP_Cox_AURKA_MYCN_all_risk,
                 exact=TRUE,
                 cens.model="cox",
                 splitMethod="none",
                 B=0,confLevel = 0.95,
                 verbose=TRUE)

plot(PredError, xlim=c(0,3000),bty='l',las=1,
     ylim = c(0, 0.3), xlab = "Days", cex.lab=3,cex.axis=3,lwd=2, smooth = T,
     legend=F)
legend('topleft',c('Reference','Subgroup',,
                   'Risk'),lwd = 2,bty='n',
       col = c(1,2, 3,4,5,6),cex=1,y.intersp=0.3)



