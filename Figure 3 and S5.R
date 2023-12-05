#### Figure 3 ####
#### Figure 3a Complex Heatmap ####
d<-as.matrix(merge_data_batch_MYCN_normal_train)
11583*0.5
mads=apply(d,1,mad) ### measured by median absolute deviation
d=d[rev(order(mads))[1:5792],] ### top 50% variable genes # 

merge_data_batch_MYCN_amp_train_train1<-data.frame(train_1_sil,train_2_sil,train_3_sil)
merge_data_batch_MYCN_amp_train_train1_heatmap<-t(merge_data_batch_MYCN_amp_train_train1[rownames(d),])
n=t(scale(merge_data_batch_MYCN_amp_train_train1_heatmap)) # 'scale'log-ratio
n[n>3]=3
n[n< -3]= -3
ac<-data.frame(c(rep('Subgroup1',1034),rep('Subgroup2',316),# Train only
                 rep('Subgroup3',370)))
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
merge_data_batch_MYCN_amp_heatmap<-t(data.frame(merge_data_batch_MYCN_amp)[rownames(d),])
n=t(scale(merge_data_batch_MYCN_amp_heatmap)) # 'scale' log-ratio
n<-n[row_order_ht,]
n[n>3]=3
n[n< -3]= -3
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

#### Supplementary Figure 5 WGCNA ####
library(WGCNA)
library(doParallel) ##
registerDoParallel(cores=4)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
merge_data_batch_MYCN_amp_train_train1_wgcna<-merge_data_batch_MYCN_amp_train_train1_heatmap
# Call the network topology analysis function
sft = pickSoftThreshold(merge_data_batch_MYCN_amp_train_train1_wgcna, 
                        powerVector = powers, verbose = 5)
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

#### WGCNA step by step ####
datExpr<-merge_data_batch_MYCN_amp_train_train1_wgcna
softPower = 6
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
pdf('WGCNA_geneTree_100.pdf',width = 100,height = 10)
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
pdf('WGCNA_MergedgeneTree_100.pdf',width = 100,height = 10)
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
pdf('WGCNA_MergedgeneTree_new_100.pdf',width = 10,height = 10)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
