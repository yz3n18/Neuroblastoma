# 02.09.2021
# generating predictions  new microarray samples; Roman Kurilov 

library(limma)
library(Biobase)

setwd("/abi/data/kurilov/koeln_agilent_custom/data2021/")

### I. collecting rawdata ###

datapath1 = "/abi/data/kurilov/koeln_agilent_custom/data2021/data_september2/"

files1 <- list.files(path = datapath1)
rawdata  = read.maimages(files=files1, path=datapath1, source="agilent",green.only=TRUE,
                             annotation= c("Row", "Col","FeatureNum", "ProbeUID", "ControlType","ProbeName", "GeneName", "SystematicName"), verbose=FALSE)
save(rawdata, file="rawdata_september2.RData") 

### II. preprocessing 

# adding raw training samples in order to normalize new samples together with training data
load("rawdata.training.RData")

normalize_roma <- function(rawdat.test){
  ncol <- ncol(rawdat.test)
  rawdat.test <- cbind(rawdat.test, rawdat.training)
  
  G.mat <- backgroundCorrect(rawdat.test, method="normexp")
  G.mat.quantnorm.test.asinh = asinh(normalizeQuantiles(G.mat$E))
  
  virtual.median = apply(G.mat.quantnorm.test.asinh, 1, median, na.rm = TRUE)
  M.quantnorm.test = G.mat.quantnorm.test.asinh-virtual.median
  #M.quantnorm.test = G.mat.quantnorm.test.asinh
  
  ctrl.spots = rawdat.test$genes$ControlType!=0 
  gene.names = rawdat.test$genes$ProbeName[!ctrl.spots]
  M.quantnorm.test = M.quantnorm.test[!ctrl.spots,]
  
  M.quantnorm.test <- as.matrix(M.quantnorm.test)
  dimnames(M.quantnorm.test) <- list(gene.names, colnames(rawdat.test$E))
  rownames(M.quantnorm.test) = sub("UKv4_", "", rownames(M.quantnorm.test))
  #colnames(M.quantnorm.test)= sub("US22502540_", "", colnames(G.mat$E))
  #colnames(M.quantnorm.test)= sub("SG14245029_", "", colnames(G.mat$E))
  
  #now let's remove training samples from normalized data
  M.quantnorm.test <- M.quantnorm.test[,1:ncol]
  return(M.quantnorm.test)
}

# now let's create an expressionset
M.quantnorm <- normalize_roma(rawdata)
eset_norm <- ExpressionSet(assayData = assayDataNew(exprs = M.quantnorm))

### III. generating predictions
library(e1071)
load("svm.model2.RData")

fs <- colnames(svm.model$SV)
svm.pred <- predict(svm.model, newdata=t(M.quantnorm[fs,]), probability = T)

write.csv(cbind(as.character(svm.pred), attr(svm.pred, "probabilities")), file="predictions.csv")



