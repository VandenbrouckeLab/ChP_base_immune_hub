## Script for downstream analysis of LPS kinetics Bulk RNA-seq
## Limma - edgeR workflow for processing, QC and DEA of HTseqCount expression table

## Load packages
library("edgeR")
library("limma")
library("ggplot2")
library("dplyr")
library("stringr")
library("openxlsx")

########################################
##### Functions
########################################

###Get DE genes
getDEgenes<-function(expMatrix, pValCutOff, logFCcutOff){
  topgenes<-expMatrix[expMatrix$adj.P.Val<pValCutOff,]
  genes_up<-topgenes[topgenes$logFC>logFCcutOff,]
  genes_down<-topgenes[topgenes$logFC< -logFCcutOff,]
  ##Sort genes on logFC
  genes_up<-genes_up[order(genes_up$logFC, decreasing=TRUE),]
  genes_down<-genes_down[order(genes_down$logFC, decreasing=TRUE),]
  genes_de_sorted<-rbind(genes_up, genes_down)
  
  return(genes_de_sorted)
}

###Get ggplot colors
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

################################################################################
######### LOAD DATA
################################################################################

getwd()
setwd("/home/clintdn/VIB/DATA/Roos/Bulk_RNA-seq_LPS/New_analysis/")
experiment<-"ChP_LPS_kinetics"

dir.create("results/")
dir.create("results/plots/")
dir.create("results/Robjects/")

##### Load raw counts
countData_all<-read.xlsx("../exp1761-RNAseqCounts.xlsx", sheet = 1)
dim(countData_all)
# 38087    23

### Preprocess data matrix (RefSeq geneID not unique!)
rownames(countData_all)<-countData_all$ID
countData<-countData_all[,4:15]
dim(countData)
# 38087    12

### Issue with non-numeric columns
for (i in 1:ncol(countData)){
  countData[,i]<-as.numeric(countData[,i])
}
countData<-as.matrix(countData)

### Make/Read meta data
colData<-read.xlsx("Metadata_Urvb_Bulk_LPS.xlsx")
rownames(colData)<-colData$SampleID
dim(colData)
# 12  3

### Reorganize and rename
cbind(colnames(countData), rownames(colData))
countData<-countData[,rownames(colData)]
cbind(colnames(countData), rownames(colData))

################################################################################
########## CREATE OBJECT
################################################################################

y <- DGEList(counts = countData)

################################################################################
########## FILTER DATA
################################################################################

##### Filter low count genes
## always work with count-per-million (CPM) instead of raw counts
## Usually a gene is required to have a count of 5-10 in a library to be considered expressed in that library
## Imagine that the lowest lib size is around 6 million reads => threshold is set on CPM>1
## But what if lib size is 20 million? Here CPM of 1 means 20 counts. Do we still use 5 counts and then the CPM cut off
## will be 0.25 or is the threshold of 5 counts a cut off for lib sizes of around 5-6 million? Then we need to put the
## cut off on 20 for lib sizes around 20 million and so use a CPM of 1.
## Threshold needs to be true in at least x samples. x is always the lowest number of replicates.
## for example: 3 samples with each 2 replicates => x set on 2
## => This ensures that a gene will be retained if it is only expressed in both replicates of a certain group

## Do filtering
yNoFilter<-y
myCpm<-cpm(y)

keep = rowSums(cpm(y)>1) >= 3 #Groups of 3!
y = y[keep,]
dim(y)
##15464
dim(yNoFilter)
##38087

##### Reset lib sizes
y$samples$lib.size = colSums(y$counts)
y$samples

################################################################################
########## NORMALISATION
################################################################################

##### Scale normalisation
yNoNormalisation<-y
y <- calcNormFactors(y)

##### Create colors for the samples
colorsUnique<-ggplotColours(4) #Number of groups
replicates<-c(3,3,3,3) #Look at group sizes

theColors<-c()
for(i in 1:length(colorsUnique)){
  theColors<-c(theColors, rep(colorsUnique[i],replicates[i]))
}

cbind(colData, theColors)

##### MDS-plot
plotMDS(y, cex=0.9, col=theColors)

png(file=paste0("results/plots/4a_MDSplot_clint_",experiment,".png"), width = 1140, height = 850)
plotMDS(y, cex=0.9, col=theColors)
dev.off()


################################################################################
########## LOG2 TRANSFORMATION
################################################################################

#### Create design matrix
colnames(colData)[3]<-"Timepoint"
TS <- paste0(colData$Condition,"_",colData$Timepoint)
TS <- factor(TS, levels=unique(TS))
design <- model.matrix(~0+TS)
colnames(design) <- levels(TS)
design

##### Do voom
png(file=paste0("results/plots/meanVariancePlot_",experiment,".png"),width = 1515, height = 1138)
v <- voom(y, design, plot = TRUE) # Transforms count data to logCPM + estimates Mean-Variance relationship (weights)
dev.off()

expTable<-v$E

### Write results
write.table(expTable, paste0("results/expTable_clint_",experiment,".txt"), sep="\t")
saveRDS(expTable, file = paste0("results/Robjects/expTable_clint_",experiment,".rds"))

## Update 2023 for upload to GEO
expTable <-readRDS( file = paste0("results/Robjects/expTable_clint_",experiment,".rds"))

cbind(colnames(countData), colnames(expTable))
all(colnames(countData) == colnames(expTable))
countData_paper<-countData
expTable_paper<-expTable

write.table(countData_paper, "Bulk_RNA_seq_raw_counts_ChP_LPS_kinetics.txt", sep="\t")
write.table(expTable_paper, "Bulk_RNA_seq_normalized_counts_ChP_LPS_kinetics.txt", sep="\t")

##### Normalised counts
countData_norm<-cpm(y, normalized.lib.sizes=TRUE, log=FALSE)

################################################################################
########## CHECK FILTERING AND NORMALISATION
################################################################################

#################### BARPLOT ####################

##### Barplot lib sizes raw counts
png(file=paste0("results/plots/1a_barplot_beforeNorm_",experiment,".png"), width = 1515, height = 1138)
par(mar = c(9,3,3,1)) #more margin: bottom, left, top, right
bp<-barplot(yNoFilter$samples$lib.size*1e-6,axisnames=FALSE,main="Barplot lib sizes of raw counts",ylab="Library size (millions)")
axis(1, labels=rownames(yNoFilter$samples), at = bp, las=2, cex.axis=0.8)
dev.off()

##### Barplot lib sizes normalised counts
png(file=paste0("results/plots/1b_barplot_afterNorm_",experiment,".png"), width = 1515, height = 1138)
par(mar = c(9,4,3,1)) #more margin: bottom, left, top, right
bp<-barplot(colSums(countData_norm)*1e-6,axisnames=FALSE,main="Barplot lib sizes of normalised counts",ylab="Library size (millions)")
axis(1, labels=colnames(countData_norm), at = bp, las=2, cex.axis=0.7)
dev.off()

#################### BOXPLOT ####################
col <- rainbow(nrow(colData))

y2<-y
y2$samples$norm.factors<-1
y2$counts[,1]<-ceiling(y2$counts[,1]*0.05)
y2$counts[,2]<-y2$counts[,2]*5

# par(mfrow=c(1,2), mar = c(15,4,3,1)) #more margin: bottom, left, top, right
png(file=paste0("results/plots/2a_boxplot_beforeNorm_",experiment,".png"), width = 1515, height = 1138)
par(mar = c(9,4,3,1)) #more margin: bottom, left, top, right
lcpm<-cpm(y2, log=T)
boxplot(lcpm, las=2, main="Unnormalised data", ylab="log-cpm", col=col)
dev.off()

png(file=paste0("results/plots/2b_boxplot_afterNorm_",experiment,".png"), width = 1515, height = 1138)
par(mar = c(9,4,3,1)) #more margin: bottom, left, top, right
y2<-calcNormFactors(y2)
lcpm<-cpm(y2, log=T)
boxplot(lcpm, las=2, main="Normalised data", ylab="log-cpm", col=col)
dev.off()

#################### DENSITY PLOT ####################
col <- topo.colors(nrow(colData))

png(file=paste0("results/plots/3a_densityPlot_",experiment,".png"), width = 1515, height = 1138)
par(mfrow=c(1,2))
### Plot log2-CPM values of each sample before filtering
theCpmNoFilter<-cpm(yNoFilter, log=TRUE)

plot(density(theCpmNoFilter[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2, main="raw data", xlab="log cpm")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0(i,". Add line for sample ",colnames(theCpmNoFilter)[i]))
  
  den <- density(theCpmNoFilter[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}

### Plot log2-CPM values of each sample after filtering (and normalisation)
theCpm<-cpm(y, log=TRUE)

plot(density(theCpm[,1]), col=col[1], lwd=2, ylim=c(0,0.25), las=2, main="filtered data", xlab="log cpm")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0(i,". Add line for sample ",colnames(theCpm)[i]))
  
  den <- density(theCpm[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}
dev.off()

par(mfrow=c(1,1))

#################### HISTOGRAM OF EXPTABLE ####################

png(file=paste0("results/plots/3b_histogramFiltering_",experiment,".png"), width = 1515, height = 1138)
par(mfrow=c(1,2))
### Histogram
hist(expTable)

### Density plot
plot(density(expTable[,1]), col=col[1], lwd=2, ylim=c(0,0.25), las=2, main="Density of expTable", xlab="log2")
abline(v=0, lty=3)
ci = 2
for (i in 2:nrow(colData)){
  print(paste0("Add line for sample ",colnames(expTable)[i]))
  
  den <- density(expTable[,i])
  lines(den$x, den$y, col=col[ci], lwd=2)
  ci = ci+1
}
dev.off()


par(mfrow=c(1,1))

################################################################################
########## PCA
################################################################################
library("rgl")

### Calculate variance
variance<-apply(expTable, 1, var)
varianceSorted<-sort(variance, decreasing=TRUE, index.return=TRUE)
### Get top 15%
numberOfGenes<-0.15*length(variance)
indexTopVariance<-varianceSorted$ix[1:numberOfGenes]
matrixPCAtmp<-expTable[indexTopVariance,]

### Prepare PCA-plot
pca<-prcomp(scale(t(matrixPCAtmp)))
matrixPCA<-cbind(pca$x[,1],pca$x[,2],pca$x[,3])
PCAcolors<-theColors

PoV <- pca$sdev^2/sum(pca$sdev^2)
summary(pca)

### Draw PCA-plot, all labels
pcaPlot<-plot3d(matrixPCA, main="",col=PCAcolors,pch=21, type="s", radius=2, legend=TRUE, xlab=paste0("pc1 (",round(PoV[1]*100,2),"%)"), 
                ylab=paste0("pc2 (",round(PoV[2]*100,2),"%)"), zlab=paste0("pc3 (",round(PoV[3]*100,2),"%)"))
text3d(x=matrixPCA[,1], y=(matrixPCA[,2]-2), z=(matrixPCA[,3]), rownames(matrixPCA) ,col=PCAcolors, cex=1)

rgl.viewpoint(0, 0)
rgl.snapshot(paste0("results/plots/4_pca_view1_",experiment,".png"))
rgl.viewpoint(35, 0)
rgl.snapshot(paste0("results/plots/4_pca_view2_",experiment,".png"))

### Save 3D
dirToSave<-paste0(getwd(),"/results/plots/")
writeWebGL(dir = dirToSave, filename = file.path(dirToSave, "4_pca.html"),
           template = system.file(file.path("WebGL", "template.html"), package = "rgl"),
           snapshot = TRUE, font = "Arial")

### Create nice plot
rgl.snapshot(paste0("results/plots/4_pca_withGrid_",experiment,".png"))


################################################################################
########## CORRELATION HEATMAP SAMPLES
################################################################################
library("RColorBrewer")
library("gplots")

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

##heatmap 1: based on distance
distsRL <- dist(t(expTable),method="euclidean")
hc <- hclust(distsRL,method="ward.D")

pdf(paste0("results/plots/4_corrSamples_distance_clint_",experiment,".pdf"))
heatmap.2(as.matrix(distsRL),
          Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col=rev(hmcol),margin=c(13, 13), cexRow=0.6, cexCol=0.6)
dev.off()

##heatmap 2: based on correlation
cm=cor(expTable)

distsRL <- as.dist(1-cor(expTable))
hc <- hclust(distsRL,method="ward.D")

pdf(paste0("results/plots/4_corrSamples_correlation_",experiment,".pdf"))
heatmap.2(cm, Rowv=as.dendrogram(hc),
          symm=TRUE, trace="none",
          col=hmcol, margin=c(13, 13), cexRow=0.6, cexCol=0.6)
dev.off()

################################################################################
########## GET DE GENES
################################################################################

head(design)
colnames(design)

#### Fit linear model on data
fit <- lmFit(expTable, design)

#### Create contrast matrix
cont.matrix <- makeContrasts(group1=LPS_1-untreated_0,
                             group2=LPS_6-untreated_0,
                             group3=LPS_24-untreated_0,
                             levels=design)
cont.matrix
fit2 = contrasts.fit(fit, cont.matrix)

#### Moderate T-test with 'shrunken' se
fit.eb <- eBayes(fit2)

#### Quick list of DE genes
summa.fit <- decideTests(fit.eb, p.value = 0.05, lfc = 1)
summary(summa.fit)

########################################
##### Extract DE genes
########################################

##### Put DE genes in list
coef<-c(1:ncol(cont.matrix))
listallgenes<-list()
listDEgenes<-list()
listDEgenes_moreStrict<-list()

for(i in 1:length(coef)){
  allGenes<-topTable(fit.eb, adjust="BH", sort.by="P",number=Inf, coef=coef[i])
  DEgenes<-getDEgenes(allGenes,0.05,1)
  DEgenes_strict<-DEgenes[DEgenes$AveExpr >0,]
  listallgenes[[i]]<-allGenes
  listDEgenes[[i]]<-DEgenes
  listDEgenes_moreStrict[[i]]<-DEgenes_strict
}


names(listallgenes)<-c("LPS_1h_vs_untreated","LPS_6h_vs_untreated","LPS_24h_vs_untreated")
names(listDEgenes)<-names(listallgenes)
names(listDEgenes_moreStrict)<-names(listDEgenes)

##### Get numbers of DE genes
lapply(listDEgenes,dim)
# 152 2629 5017

########################################
##### Write results
########################################

### Add geneSymbol in column (for the export)
listDEgenes<-lapply(listDEgenes,function(x){dplyr::mutate(x,'gene'=rownames(x))})
listallgenes<-lapply(listallgenes,function(x){dplyr::mutate(x,'gene'=rownames(x))})

### Write results
library('openxlsx')
write.xlsx(listDEgenes, file = paste0("results/summary_DEgenes_clint_",experiment,".xlsx"))
write.xlsx(listallgenes, file = paste0("results/summary_allgenes_clint_",experiment,".xlsx"))

### Save results
saveRDS(listDEgenes, file=paste0("results/Robjects/listDEgenes_clint_",experiment,".rds"))
saveRDS(listallgenes, file=paste0("results/Robjects/listallgenes_clint_",experiment,".rds"))
