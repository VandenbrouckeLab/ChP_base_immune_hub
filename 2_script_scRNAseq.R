## Script for processing scRNA-seq data
## Analysis from 10X CellRanger data matrices 
## Contains pre-processing, QC, Seurat workflow and all custom visualizations from the paper

library('Seurat')
library('dplyr')
library('gridExtra')
library('scater')
library('DisneyTools')
library('RColorBrewer')
library('ggrepel')
library("ggthemes")
library("scales")
library("ggpubr")
library('gplots')
library('openxlsx')
library('future')
plan("multiprocess", workers = 8)

library('DoubletFinder')
library('fields')
library('modes')
library("dplyr")
library("stringr")

library('clusterProfiler')
library('org.Mm.eg.db')

################################################################################
########## GENERAL
################################################################################

########################################
##### Getwd
########################################

setwd("~/VIB/DATA/Roos/Daan 1/")

sampleName<-"RVD_aggr17"
sampleFolder<-paste0(sampleName,"/")

########################################
##### Some variables
########################################

### Read from the file "aggregation.csv"
aggrFile<-read.csv(file=paste0(sampleFolder,"aggregation.csv"), stringsAsFactors = F)
if(ncol(aggrFile)==3){
  listLabels<-as.list(paste0(aggrFile[,1],"-",aggrFile[,3]))
}else{
  listLabels<-as.list(aggrFile[,1])
}
listLabels

### General variables
diagnostics<-list()

########################################
##### Functions
########################################
source('~/VIB/DATA/Roos/Daan 1/script_functions.R')
source('~/VIB/DATA/Roos/Daan 1/script_featurePlots.R')

################################################################################
########## LOAD DATA
################################################################################
rawDataSparse <- Read10X(paste0(sampleFolder,"filtered_feature_bc_matrix/"))
dim(rawDataSparse)
# [1] 27998 23782

########################################
##### Characteristics of data
########################################

nrZeros<-sum(rawDataSparse==0)/(nrow(rawDataSparse)*ncol(rawDataSparse))*100
nrZeros
## 94.56197

##### In each cell: how many genes are expressed (count > 0) #####
cellCounts<-apply(rawDataSparse,2,function (x){sum(x>0)})
length(cellCounts[cellCounts<200])
## 163

##### For each gene: in how many cells is it expressed? (count > 0) #####
geneCounts<-apply(rawDataSparse,1,function (x){sum(x>0)})
length(geneCounts[geneCounts<3])
## 9985

##### Add to diagnostics #####
diagnostics[['dimrawDataSparse']]<-paste0(nrow(rawDataSparse)," genes - ",ncol(rawDataSparse)," cells")
diagnostics[['nrGenes']]<-nrow(rawDataSparse)
diagnostics[['nrCells']]<-ncol(rawDataSparse)
diagnostics[['zeroInflation']]<-nrZeros

diagnostics[['minGenesPerCell']]<-min(cellCounts)
diagnostics[['maxGenesPerCell']]<-max(cellCounts)
diagnostics[['meanGenesPerCell']]<-mean(cellCounts)
diagnostics[['medianGenesPerCell']]<-median(cellCounts)
diagnostics[['cellsLess200genes']]<-length(cellCounts[cellCounts<200])
diagnostics[['genesNotExpressed']]<-length(geneCounts[geneCounts==0])
diagnostics[['genesLess3cells']]<-length(geneCounts[geneCounts<3])

### Remove some variables
rm(geneCounts)
rm(cellCounts)

################################################################################
########## QC: CELLS
################################################################################

########################################
########## Calculate metrics
########################################

##### Create object #####
library("scater")
sce<-SingleCellExperiment(list(counts=rawDataSparse))
dim(sce)
diagnostics[['dimSce']]<-paste0(nrow(sce)," genes - ",ncol(sce)," cells")

##### Get spike inns #####
is.spike <- grepl("^ercc-", rownames(sce))
sum(is.spike)
##0

##### Get mitochondrial genes #####
is.mito <- grepl("^mt-", rownames(sce))
sum(is.mito)
##13
rownames(sce)[is.mito]

##### Calculate QC metrics #####
### => pData(sce) is created
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito))
dim(colData(sce))

##### Create metaData matrix (used for downstream analysis) #####
metaData<-data.frame("staticNr"=colnames(rawDataSparse),"orig.ident"=listLabels[[1]],
                     "nGene"=sce$total_features_by_counts,"nUMI"=sce$total_counts,"percent.mito"=sce$pct_counts_Mt, 
                     stringsAsFactors = F)
rownames(metaData)<-metaData$staticNr
metaData$staticNr<-1

for(i in 2:length(listLabels)){
  toSearch<-paste0("-",i)
  metaData[grep(toSearch,rownames(metaData)), which(colnames(metaData)=="orig.ident")]<-listLabels[[i]]
}

table(metaData$orig.ident)

##### Add to diagnostics #####
diagnostics[['splitSamples']]<-paste0(table(metaData$orig.ident)," cells of sample ",rownames(table(metaData$orig.ident)))

########################################
########## Get outliers
########################################
nmad_low_feature<-4
nmad_high_feature<-4

nmad_low_UMI<-4
nmad_high_UMI<-4

nmad_high_mito<-5

##### Aim: remove cells with low library sizes, low numbers of expressed features and with high mitochondrial proportions
##same as nGene in Seurat pipeline
feature.drop.low <- isOutlier(sce$total_features_by_counts, nmads=nmad_low_feature, type="lower", log=TRUE)
sum(feature.drop.low)

feature.drop.high <- isOutlier(sce$total_features_by_counts, nmads=nmad_high_feature, type="higher", log=TRUE)
sum(feature.drop.high)

feature.drop<-as.logical(feature.drop.low + feature.drop.high)
sum(feature.drop)

##same as UMI in Seurat pipeline
libsize.drop.low <- isOutlier(sce$total_counts, nmads=nmad_low_UMI, type="lower", log=TRUE)
sum(libsize.drop.low)

libsize.drop.high <- isOutlier(sce$total_counts, nmads=nmad_high_UMI, type="higher", log=TRUE)
sum(libsize.drop.high)

libsize.drop<-as.logical(libsize.drop.low+libsize.drop.high)
sum(libsize.drop)

##% mitochondrial genes
mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=nmad_high_mito, type="higher")
sum(mito.drop) #1061

##### add to metaData matrix #####
metaData$nGene.drop=feature.drop
metaData$nUMI.drop=libsize.drop
metaData$mito.drop=mito.drop
metaData$final.drop=feature.drop | libsize.drop | mito.drop

########################################
########## Create histogram + barplot
########################################
palette(c("#00BFC4","#F8766D","#7CAE00","#C77CFF"))

toPlot<-metaData
savePlots<-TRUE

##nGene
if(savePlots==TRUE){png(file=paste0(sampleFolder,"results/QC/1a_nGene.png"), width=850)}
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nGene),]
hist(tmp$nGene, breaks=30)
theColors<-as.factor(tmp$nGene.drop)
barplot(tmp$nGene, col=theColors, border=theColors)
if(savePlots==TRUE){dev.off()}

##nUMI
if(savePlots==TRUE){png(file=paste0(sampleFolder,"results/QC/1b_nUMI.png"), width=850)}
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$nUMI),]
hist(tmp$nUMI, breaks=30)
theColors<-as.factor(tmp$nUMI.drop)
barplot(tmp$nUMI, col=theColors, border=theColors)
if(savePlots==TRUE){dev.off()}

##percent.mito
if(savePlots==TRUE){png(file=paste0(sampleFolder,"results/QC/1c_percMito.png"), width=850)}
par(mfrow=c(1,2))
tmp<-toPlot[order(toPlot$percent.mito),]
hist(tmp$percent.mito, breaks=30)
theColors<-as.factor(tmp$mito.drop)
barplot(tmp$percent.mito, col=theColors, border=theColors)
if(savePlots==TRUE){dev.off()}

########################################
########## Create violinPlots
########################################

### Before filtering
toPlot<-metaData
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/2a_beforeFiltering.png"), colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))
drawVlnPlot_split(toPlot, fileName = paste0(sampleFolder,"results/QC/2a_beforeFiltering_splitted.png"), colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))

drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/2a_beforeFiltering_nGene.png"), 
            colsToColor = c('nGene.drop','nGene.drop','nGene.drop'))
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/2a_beforeFiltering_nUMI.png"), 
            colsToColor = c('nUMI.drop','nUMI.drop','nUMI.drop'))
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/2a_beforeFiltering_mito.png"), 
            colsToColor = c('mito.drop','mito.drop','mito.drop'))

### After filtering
toPlot<-metaData[! metaData$final.drop,]
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/2b_afterFiltering.png"), colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))
drawVlnPlot_split(toPlot, fileName = paste0(sampleFolder,"results/QC/2b_afterFiltering_splitted.png"), colsToColor = c('nGene.drop','nUMI.drop','mito.drop'))

########################################
########## Remove outliers
########################################

sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
dim(sce)

### Number of cells removed
nrow(metaData)-ncol(sce) #1071

########################################
########## Create PCA 
########################################
library('mvoutlier')

varsToUse <- c("pct_counts_in_top_100_features",
               "total_features_by_counts", "pct_counts_feature_control",
               "total_features_by_counts_feature_control", "log10_total_counts_endogenous",
               "log10_total_counts_feature_control")
setdiff(varsToUse, colnames(colData(sce)))
exprs_to_plot <- scale(colData(sce)[,varsToUse], scale = T)
x.mad = apply(exprs_to_plot, 2, mad)
x.mad[x.mad==0]
varsToUse<-setdiff(varsToUse, names(x.mad[x.mad==0]))

##### Detect bad cells #####
sceNew<-runColDataPCA(sce,outliers = T, variables = varsToUse)
table(sceNew$outlier)
#1251 cells dropped as outliers

outs<-colnames(sceNew)[sceNew$outlier]

### Add to metaData
metaData$pca.drop<-metaData$final.drop
metaData[outs,which(colnames(metaData)=="pca.drop")]<-TRUE

##### Color bad cells on PCA plot #####
colorDF<-as.data.frame(cbind(colnames(sceNew),"1"), stringsAsFactors=F)
rownames(colorDF)<-colorDF[,1]
colorDF[outs,2]<-"2"
colorDF[,2]<-as.factor(colorDF[,2])
tmp<-colorDF[,2,drop=F]

png(file=paste0(sampleFolder,"results/QC/3a_pca.png"),  width = 850, height = 642)
plotReducedDim(sceNew, dimred = "PCA_coldata", colour_by='outlier',shape_by='outlier') + labs(title="PCA with outliers colored")
dev.off()

#### Add to metaData table ####
pca.drop<-metaData[colnames(sce),"pca.drop"]
sum(pca.drop)

##### Create violinplots ####
##Before
toPlot<-metaData[! metaData$final.drop,]
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/3b_beforePcaFiltering_2.png"), colsToColor = c(rep('pca.drop',3)))

##After
toPlot<-metaData[! metaData$pca.drop,]
drawVlnPlot(toPlot, fileName = paste0(sampleFolder,"results/QC/3c_afterPcaFiltering.png"), colsToColor = c(rep('pca.drop',3)))

##### Remove outlier cells #####
sce <- sce[,!(pca.drop)]
dim(sce)
diagnostics[['dimAfterPCA']]<-paste0(nrow(sce)," genes - ",ncol(sce)," cells")

##### Add to diagnostics #####
diagnostics[['pcaRemove']]<-0
diagnostics[['pcaRemove']]<-sum(pca.drop)
diagnostics[['totalRemove']]<-nrow(metaData)-ncol(sce)

### Remove some variables
rm(sceNew)

################################################################################
########## FINALIZE QC
################################################################################

dim(sce)
saveRDS(sce, file=paste0(sampleFolder,"Robjects/sce.rds"))

rawDataSparseFiltered<-rawDataSparse[rownames(sce),colnames(sce)]
dim(rawDataSparseFiltered)
# 27998 21460
diagnostics[['dimBeforeSeuratObj']]<-paste0(nrow(rawDataSparseFiltered)," genes - ",ncol(rawDataSparseFiltered)," cells")

### Remove some variables
rm(sce)
rm(rawDataSparse)

################################################################################
########## CREATE SEURAT OBJECT
################################################################################

##### Create object #####
seuratObj <- CreateSeuratObject(counts = rawDataSparseFiltered, project = "seuratObj", min.cells = 3, min.features = 200)

### Explore object
dim(seuratObj)
#17815 21460

##### Add to diagnostics #####
diagnostics[['dimAfterSeuratObj']]<-paste0(nrow(seuratObj)," genes - ",ncol(seuratObj)," cells")

################################################################################
########## FILTER DATA
################################################################################

seuratObj[["percent.mito"]] <- PercentageFeatureSet(object = seuratObj, pattern = "^mt-")
head(seuratObj@meta.data)

png(file=paste0(sampleFolder,"results/QC/4_vlnPlotSeurat.png"), width = 850, height = 642)
VlnPlot(object = seuratObj, features = c("nFeature_RNA", "nCount_RNA","percent.mito"), ncol = 3)
dev.off()

##### Add orig ident
metaDataTable<-seuratObj@meta.data
metaDataTable$orig.ident<-as.character(metaDataTable$orig.ident)
for(i in 1:length(listLabels)){
  toSearch<-paste0('-',i)
  metaDataTable[grep(toSearch,rownames(metaDataTable)), which(colnames(metaDataTable)=="orig.ident")]<-listLabels[[i]]
}
seuratObj@meta.data<-metaDataTable
head(metaDataTable)
table(metaDataTable$orig.ident)

##### Add to diagnostics #####
diagnostics[['splitSamplesAfterFiltering']]<-paste0(table(metaDataTable$orig.ident)," cells of sample ",rownames(table(metaDataTable$orig.ident)))

################################################################################
########## NORMALIZE
################################################################################
seuratObj <- NormalizeData(object = seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)

##### Check per group #####
metaDataTable<-seuratObj@meta.data
metaDataTable$nUMI<-colSums(as.matrix(seuratObj[['RNA']]@data))
metaDataTable$nGene<-apply(as.matrix(seuratObj[['RNA']]@data),2,function(x){sum(x>0)})

drawVlnPlotSeurat_split(metaDataTable, paste0(sampleFolder,"results/QC/5_afterNorm_splitted.png"))

################################################################################
########## GET HVG
################################################################################

seuratObj <- FindVariableFeatures(object = seuratObj, selection.method = "vst", nfeatures = 2000)
length(VariableFeatures(seuratObj))

### Get more info about HVGs (mean, dispersion and dispersion scaled)
head(HVFInfo(seuratObj))

### Plot variable features with and without labels
top10 <- head(x = VariableFeatures(object = seuratObj), 10)
plot1 <- VariableFeaturePlot(object = seuratObj)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

png(file=paste0(sampleFolder,"results/QC/6_hvg.png"), width = 850, height = 642)
CombinePlots(plots = list(plot1, plot2))
dev.off()

##### Add to diagnostics #####
diagnostics[['varGenes']]<-length(VariableFeatures(seuratObj))

################################################################################
########## SCALE DATA
################################################################################
# Apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
# Scaling data => Shifts the expression of each gene, so that the mean expression across cells is 0
# Scaling data => Scales the expression of each gene, so that the variance across cells is 1 (so that highly-expressed genes do not dominate)

seuratObj <- ScaleData(object = seuratObj, features = rownames(seuratObj))

##### Check per group #####
head(metaDataTable)
metaDataTable$nUMI<-colSums(as.matrix(seuratObj[['RNA']]@scale.data))
metaDataTable$nGene<-apply(as.matrix(seuratObj[['RNA']]@scale.data),2,function(x){sum(x>0)})

drawVlnPlotSeurat_split(metaDataTable, paste0(sampleFolder,"results/QC/7_afterScale_splitted.png"))

################################################################################
########## PCA
################################################################################
seuratObj <- RunPCA(object = seuratObj, features = VariableFeatures(seuratObj), npcs = 50, ndims.print = 1:5, nfeatures.print = 10)

seuratObj[['pca']]@cell.embeddings[1:5,1:5]

########################################
########## PCA PLOT
########################################
pdf(file=paste0(sampleFolder,"results/QC/8a_PCA.pdf"), width = 10)
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,3))
dev.off()

########################################
########## PCA PLOT 3D
########################################
library("rgl")

expTable<-seuratObj[['RNA']]@data
matrixPCAtmp<-expTable[VariableFeatures(seuratObj),]

### Prepare PCA-plot
pca<-seuratObj[['pca']]@cell.embeddings
matrixPCA<-cbind(pca[,1],pca[,2],pca[,3])

PoV <- seuratObj[['pca']]@stdev^2/sum(seuratObj[['pca']]@stdev^2)
summary(pca[,c(1,2,3)])

### Draw PCA-plot, all labels
pcaPlot<-plot3d(matrixPCA, main="",pch=2, type="s",radius=0.5, legend=TRUE, xlab=paste0("pc1 (",round(PoV[1]*100,2),"%)"), 
                ylab=paste0("pc2 (",round(PoV[2]*100,2),"%)"), zlab=paste0("pc3 (",round(PoV[3]*100,2),"%)"))

rgl.viewpoint(0, 5)
rgl.snapshot(paste0(sampleFolder,"results/QC/8b_PCA_view1.png"))
rgl.viewpoint(35, 5)
rgl.snapshot(paste0(sampleFolder,"results/QC/8b_PCA_view2.png"))

########################################
########## HEATMAP OF PCs
########################################

### Create heatmap of PC 1-40
pdf(file=paste0(sampleFolder,"results/QC/9a_selectPC.pdf"))
PCHeatmap(seuratObj, dims = 1:12, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 13:24, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 25:36, cells = 500, balanced = TRUE)
PCHeatmap(seuratObj, dims = 37:40, cells = 500, balanced = TRUE)
dev.off()

################################################################################
########## DETERMINE STATISTICALLY SIGNIFICANT PCs
################################################################################

### Create PCElbowplot
png(file=paste0(sampleFolder,"results/QC/9b_selectPC.png"), width = 850, height = 642)
ElbowPlot(object = seuratObj, ndims = 40)
dev.off()

################################################################################
########## AUTOMATICALLY CLUSTER THE CELLS (35 dims res 1.0)
################################################################################

## Parameters to try 
dimsToTry<-c(25,30,35)
resToUse<-0.8

### Final parameters chosen
dimsToTry<-c(35)
resToUse<-1.0
diagnostics[['dimsPC']]<-dimsToTry
diagnostics[['res']]<-resToUse


for(maxPCs in dimsToTry){
  dimsToUse<-1:maxPCs
  print(paste0("Working on 1:",maxPCs))
  
  ##### Find clusters
  seuratObj <- FindNeighbors(object = seuratObj, dims = dimsToUse)
  seuratObj <- FindClusters(object = seuratObj, resolution = resToUse)
  
  ##### Create tSNE plot
  seuratObj <- RunTSNE(object = seuratObj, dims = dimsToUse)
  tsnePlot<-DimPlot(seuratObj, reduction = "tsne", label=T, label.size = 8, pt.size = 2)
  tsnePlotSplit<-DimPlot(seuratObj, reduction = "tsne", label=F, group.by="orig.ident", pt.size = 2)
  
  ggsave(grid.arrange(tsnePlot, tsnePlotSplit, ncol=2),
         file=paste0(sampleFolder,"results/QC/10a_tSNE_final_",min(dimsToUse),"-",max(dimsToUse),"-", resToUse,"_sliced.png"), width = 20)
  
  ##### Create UMAP plot
  seuratObj <- RunUMAP(seuratObj, dims = dimsToUse, n.neighbors = 30) 
  umapPlot<-DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)
  umapPlotSplit<-DimPlot(seuratObj, reduction = "umap", label = F, group.by="orig.ident")
  
  ggsave(grid.arrange(umapPlot, umapPlotSplit, ncol=2),
         file=paste0(sampleFolder,"results/QC/10b_UMAP_final_",min(dimsToUse),"-",max(dimsToUse),"-", resToUse,"_sliced.png"), width = 20)
  
}

################################################################################
########## DOUBLET FINDER (multiple samples)
################################################################################

minPCT = 1
maxPCT = 10
pN = 0.25
PCs=1:40
DFPredictions <- c()
nSamples <- str_split_fixed(rownames(seuratObj@meta.data), "-", 2)[,2] %>% unique() %>% length()
seuratObjDoubletFinder<-seuratObj

##### Run DoubletFinder #####
for(i in 1:nSamples) {
  cellsToUse <- grep(paste0("-", i, "$"), rownames(seuratObjDoubletFinder@meta.data), value = T)
  object_sub <- SubsetData(seuratObjDoubletFinder, cells = cellsToUse)
  nDoubletsMin <- (length(colnames(seuratObjDoubletFinder)) * (minPCT/100) ) %>% round()
  nDoubletsMax <- (length(colnames(seuratObjDoubletFinder)) * (maxPCT/100) ) %>% round()
  
  findpK <- paramSweep_v3(object_sub, PCs = PCs, sct=FALSE, num.cores = 4) %>%
    summarizeSweep() %>%
    find.pK()
  
  maxScore <- findpK %>% pull('BCmetric') %>% which.max()
  pKValue <- findpK[maxScore, 'pK'] %>% as.character() %>% as.numeric()
  
  object_sub <- doubletFinder_v3(seu = object_sub, pN = pN, pK = pKValue, nExp = nDoubletsMin, reuse.pANN = F, PCs = PCs, sct=FALSE)
  object_sub <- doubletFinder_v3(seu = object_sub, pN = pN, pK = pKValue, nExp = nDoubletsMax, reuse.pANN = paste0("pANN_", pN, "_", pKValue, "_", nDoubletsMin), PCs = PCs, sct=FALSE)
  
  object_sub@meta.data$DFPrediction <- object_sub@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_", nDoubletsMax)]
  object_sub@meta.data$DFPrediction[ object_sub@meta.data[, paste0("DF.classifications_", pN, "_", pKValue, "_", nDoubletsMin)] == 'Doublet'] <- "High Confidence"
  object_sub@meta.data$DFPrediction <- gsub('Doublet', "Low Confidence", object_sub@meta.data$DFPrediction)
  
  newData <- object_sub@meta.data$DFPrediction
  names(newData) <- rownames(object_sub@meta.data)
  
  DFPredictions <- c(DFPredictions, newData)
}

seuratObjDoubletFinder@meta.data$DFPrediction <- DFPredictions[rownames(seuratObj@meta.data)]

##### Visualize results #####
p<-DimPlot(seuratObjDoubletFinder, reduction.use = 'umap', group.by = 'DFPrediction', 
           cols.use = c("red", "yellow", "#C9C9C9"))
ggsave(p, file=paste0(sampleFolder,"results/QC/17.doubletFinder.png"))

################################################################################
########## CHECK DE GENES 
################################################################################
dir.create(paste0(sampleFolder,"results/QC/Feature_plots"))

##### Epithelial marker -> 1,2,5 + 3,4,12 <-> 19,22 Xist?? <-> 23 (doublets!!)
F1 <- FeaturePlot(object = seuratObj, features = c("Otx2", "Ttr"), cols = c("grey", "blue"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Epithelial_Cells.png"), height = 10, width = 20, dpi = "retina")

##### Endothelial marker -> 8,6 <-> little blob above 7 (doublets) SLICE!! + 30 (doublets)
F1<-FeaturePlot(object = seuratObj, features = c("Pecam1", "Flt1", "Plvap"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Endothelial_Cells.png"), height = 10, width = 20, dpi = "retina")

##### Vascular associated marker -> 17 (a split cluster) -> SLICE!!
F1<-FeaturePlot(object = seuratObj, features = c("Pdgfrb", "Mylk", "Myh11", "Tagln"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Vasc_Assoc_Cells.png"), height = 10, width = 20, dpi = "retina")

##### Fibroblast marker -> 9+10 <-> 31 (doublets)
F1<-FeaturePlot(object = seuratObj, features = c("Dcn", "Col1a1"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Fibroblasts.png"), height = 10, width = 20, dpi = "retina")

##### Macrophage marker -> 7+11+0+15+16 <-> 20 (doublets)
F1<-FeaturePlot(object = seuratObj, features = c("Adgre1", "Csf1r","Fcgr1"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Macrophages.png"), height = 10, width = 20, dpi = "retina")

##### Microglia marker -> 26 
F1 <- FeaturePlot(object = seuratObj, features = "P2ry12", cols = c("grey", "blue"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_Microglia like.png"), height = 10, width = 10, dpi = "retina")


##### NK cell marker -> left side 24 + 14 -> SLICE right side 14
F1<- FeaturePlot(object = seuratObj, features = c("Klrb1c", "Gzmb"), cols = c("grey", "blue"), 
                 reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_NK_Cells.png"), height = 10, width = 20, dpi = "retina")

##### cDC marker -> 27 Xcr1+ + 28 Cd209a+ Top half of 18 (Ccr7+)
F1<-FeaturePlot(object = seuratObj, features = c("Cd209a", "Ccr7","Xcr1"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_cDC.png"), height = 20, width = 20, dpi = 250)


##### Neutrophil marker -> 25
F1<-FeaturePlot(object = seuratObj, features = c("S100a8", "Ngp", "Retnlg"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_neutrophils.png"), height = 20, width = 20, dpi = 250)

##### Monocyte marker -> 13 + right part of 24?? + 25 shows activity??
F1<-FeaturePlot(object = seuratObj, features = c("Plac8", "Ly6c2", "Ccr2"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_monocytes.png"), height = 20, width = 20, dpi = 250)

##### Mast cell marker -> 25 + 13
F1<-FeaturePlot(object = seuratObj, features = "Mcemp1", cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_mast_cells.png"), height = 10, width = 10, dpi = 250)


##### Neuronal cell marker -> bottom left of 21 (except furthest bottom) -> SLICE twice!!
FeaturePlot(object = seuratObj, features = "Tubb3", cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

##### Mitotic cell marker -> 21 (except bottom left blob)
F1<-FeaturePlot(object = seuratObj, features = c("Birc5", "Mki67"), cols = c("grey", "blue"), 
                reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)
ggsave(F1, file=paste0(sampleFolder,"results/QC/Feature_plots/Feature_plot_mitotic_cells.png"), height = 10, width = 20, dpi = 250)

##### Glial marker -> ??
FeaturePlot(object = seuratObj, features = c("Slc1a3", "Fabp7"), cols = c("grey", "blue"), 
            reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', pt.size = 1.5)

########################################
##### First marker analysis
########################################

### Find markers for every cluster compared to all remaining cells, report only the positive ones
allMarkers <- FindAllMarkers(seuratObj, min.pct = 0.10, min.diff.pct=0.25, logfc.threshold = 0.30, return.thresh = 0.01, only.pos = TRUE)
table(allMarkers$cluster)
saveRDS(allMarkers, file=paste0(sampleFolder,"Robjects/allMarkers_sliced_",sampleName,".rds"))

### Add to diagnostics
diagnostics<-readRDS(file=paste0(sampleFolder,"Robjects/diagnostics_sliced_",sampleName,".rds"))
diagnostics[['markersPerCluster']]<-paste0(table(allMarkers$cluster)," markers for cluster ",rownames(table(allMarkers$cluster)))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_sliced_",sampleName,".rds"))

### Create list with markers
totalNrClusters<-max(as.numeric(names(table(allMarkers$cluster))))
totalNrClustersPlusOne<-totalNrClusters+1
markersList<-list()

for(i in 1:totalNrClustersPlusOne){
  clusterNr<-i-1
  
  tmp<-allMarkers[allMarkers$cluster==clusterNr,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  markersList[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}
names(markersList)<-paste0("cluster",0:totalNrClusters)

### Write to Excel
library('openxlsx')
write.xlsx(markersList, file =paste0(sampleFolder, "results/allClusters_sliced_",sampleName,".xlsx"))

################################################################################
########## EXTRA MANUAL CLUSTERING
################################################################################

##### Manually split clusters based on dimreduc and marker expression!
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

# Doublets 7
tsneSlice<-tsneTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., tSNE_1 > -33.5, tSNE_2 < -5)
wantedCells<-intersect(tsneSlice$cell, WhichCells(seuratObj, idents = 7))
colorSomeCells(clusterMatrix, tsneTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 32)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

# Vascular Associated Cells 17
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_2 < 5)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 17))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 33)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

# T cells 14
umapSlice1<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 > 2, UMAP_2 >2.8)
umapSlice2<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 > 3, UMAP_2 <2.8)
umapSlice<-rbind(umapSlice1,umapSlice2)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 14))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 34)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

# Neuron 21
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 < 2.8, UMAP_1 > 1.5)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 21))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 35)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

# Unknown 21
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 < 1.5)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 21))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 36)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

# FB2 9
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 < 5.7, UMAP_2 > -11.6)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 9))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 37)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

# FB2 10
umapSlice<-umapTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., UMAP_1 > -6.35, UMAP_1 < -4, UMAP_2 > 10.98)
wantedCells<-intersect(umapSlice$cell, WhichCells(seuratObj, idents = 10))
colorSomeCells(clusterMatrix, umapTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 38)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

# Doublets 6
tsneSlice<-tsneTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., tSNE_1 < -8, tSNE_2 > -31.8)
wantedCells<-intersect(tsneSlice$cell, WhichCells(seuratObj, idents = 6))
colorSomeCells(clusterMatrix, tsneTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 39)
DimPlot(seuratObj, reduction = "umap", label = T, label.size = 8)

# Doublets 11
tsneSlice<-tsneTable %>% dplyr::mutate('cell'=rownames(.)) %>% dplyr::filter(., tSNE_1 > -32.5, tSNE_2 > 5)
wantedCells<-intersect(tsneSlice$cell, WhichCells(seuratObj, idents = 11))
wantedCells2<-intersect(tsneSlice$cell, WhichCells(seuratObj, idents = 7))
wantedCells<-c(wantedCells,wantedCells2)
colorSomeCells(clusterMatrix, tsneTable, wantedCells)

seuratObj<-SetIdent(object = seuratObj, cells = wantedCells, value = 40)
DimPlot(seuratObj, reduction = "tsne", label = T, label.size = 8)

# Increase cluster 17 
seuratObj2 <- CellSelector(U1, object=seuratObj, ident=17)
seuratObj<-seuratObj2
U1 <- DimPlot(seuratObj, reduction = "umap", label = T, label.size = 4)

##### Save object
saveRDS(seuratObj, file=paste0(sampleFolder,"Robjects/seuratObj_sliced_",sampleName,".rds"))
saveRDS(diagnostics, file=paste0(sampleFolder,"Robjects/diagnostics_sliced_",sampleName,".rds"))

##### Read object
seuratObj <- readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_sliced_",sampleName,".rds"))
diagnostics <- readRDS(file=paste0(sampleFolder,"Robjects/diagnostics_sliced_",sampleName,".rds"))

################################################################################
########## PLOTS
################################################################################
clusterMatrix<-seuratObj@meta.data
tsneTable<-as.data.frame(seuratObj[['tsne']]@cell.embeddings, stringsAsFactors = F)
umapTable<-as.data.frame(seuratObj[['umap']]@cell.embeddings, stringsAsFactors = F)

########## UMI plot ##########
p1<-drawUMI_mitoPlot(tsneTable, 'tsne', clusterMatrix, 'nCount_RNA',"UMI")
p2<-drawUMI_mitoPlot(umapTable, 'umap', clusterMatrix, 'nCount_RNA',"UMI")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results/QC/11a_UMI.png"), width = 20)

########## mito.genes plot ##########
p1<-drawUMI_mitoPlot(tsneTable, 'tsne', clusterMatrix, 'percent.mito',"mito")
p2<-drawUMI_mitoPlot(umapTable, 'umap', clusterMatrix, 'percent.mito',"mito")

ggsave(grid.arrange(p1, p2, ncol=2), file=paste0(sampleFolder,"results/QC/11b_percMito.png"), width = 20)

########## PCA plot ##########
pdf(file=paste0(sampleFolder,"results/QC/13a_PCA.pdf"), width=10)
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,2))
DimPlot(object = seuratObj, reduction = "pca", dims = c(2,3))
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,3))
dev.off()

pdf(file=paste0(sampleFolder,"results/QC/13b_PCA_split.pdf"), width=10)
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,2), group.by = "orig.ident")
DimPlot(object = seuratObj, reduction = "pca", dims = c(2,3), group.by = "orig.ident")
DimPlot(object = seuratObj, reduction = "pca", dims = c(1,3), group.by = "orig.ident")
dev.off()

### Create initial annotated UMAP
seuratObj@meta.data$annotated_clusters <- seuratObj@active.ident
seuratObj@meta.data$annotated_clusters<- factor(seuratObj@meta.data$annotated_clusters,levels(seuratObj@meta.data$annotated_clusters)[c(11:27,1,28:41,10,9,8,7,6,5,4,3,2)]) #reorder levels
levels(seuratObj@meta.data$annotated_clusters) <- c("Macrophages","Epithelial Cells","Epithelial Cells","Epithelial Cells",
                                                    "Epithelial Cells",'Epithelial Cells','Endothelial Cells','Macrophages',
                                                    "Endothelial Cells",'Stromal Fibroblasts',"Stromal Fibroblasts","Macrophages",'Epithelial Cells','Monocytes','NK Cells','Macrophages',
                                                    "Macrophages","Vascular Associated Cells","Dendritic Cells",'Fras1+ Epithelial Cells','Doublets','Mitotic Cells','Fras1+ Epithelial Cells',
                                                    "Doublets_2","NK Cells","Neutrophils","Epiplexus Cells",
                                                    'Dendritic Cells',"Dendritic Cells",'Fibroblast-like Cells',"Doublets_3","Doublets_4",
                                                    "Doublets_5","Vascular Associated Cells","T Cells","Neuronal and Glial Cells","Neuronal and Glial Cells",
                                                    'Stalk Fibroblasts',"Stalk Fibroblasts","Doublets_6","Doublets_7") 
                                                    #DCs together, neuronal and glial together, rename FB, rename Xist+!!!!
U1 <- DimPlot(seuratObj, reduction = "umap", label = T, repel = T, label.size = 4, group.by="annotated_clusters")

### Create new clusters: split on treatment
seuratObjNew<-seuratObj
Idents(seuratObjNew)<-seuratObjNew@meta.data$annotated_clusters
seuratObjNew@meta.data$newClustersTmp<-Idents(seuratObjNew)
seuratObjNew@meta.data$Treatment<-seuratObjNew@meta.data$orig.ident
seuratObjNew@meta.data[which(seuratObjNew@meta.data$Treatment == "LpsNegFour" | seuratObjNew@meta.data$Treatment == "LpsNegLat"),"Treatment"]<-"LpsNeg"
seuratObjNew@meta.data[which(seuratObjNew@meta.data$Treatment == "LpsPosFour" | seuratObjNew@meta.data$Treatment == "LpsPosLat"),"Treatment"]<-"LpsPos"
seuratObjNew@meta.data$newClusters<-paste0(seuratObjNew@meta.data$newClustersTmp,"_",seuratObjNew@meta.data$Treatment)
head(seuratObjNew@meta.data)

### Use the new clusters
Idents(seuratObjNew)<-seuratObjNew@meta.data$newClusters

# Remove the bad clusters!!
seuratObjNew<-subset(seuratObjNew, idents =c("Epithelial Cells_LpsPos","Epithelial Cells_LpsNeg","Fras1+ Epithelial Cells_LpsPos","Fras1+ Epithelial Cells_LpsNeg",
                                             "Endothelial Cells_LpsPos","Endothelial Cells_LpsNeg","Vascular Associated Cells_LpsPos","Vascular Associated Cells_LpsNeg",
                                             "Stromal Fibroblasts_LpsPos","Stromal Fibroblasts_LpsNeg","Stalk Fibroblasts_LpsPos","Stalk Fibroblasts_LpsNeg",
                                             "Fibroblast-like Cells_LpsPos","Fibroblast-like Cells_LpsNeg","Macrophages_LpsPos","Macrophages_LpsNeg","Epiplexus Cells_LpsNeg", 
                                             "NK Cells_LpsPos","NK Cells_LpsNeg","Dendritic Cells_LpsPos","Dendritic Cells_LpsNeg",
                                             "T Cells_LpsPos","T Cells_LpsNeg","Monocytes_LpsPos","Monocytes_LpsNeg","Neutrophils_LpsPos","Mitotic Cells_LpsPos","Mitotic Cells_LpsNeg",
                                             "Neuronal and Glial Cells_LpsPos","Neuronal and Glial Cells_LpsNeg")   )

seuratObjNew@meta.data$newClusters<-factor(seuratObjNew@meta.data$newClusters, levels=c("Epithelial Cells_LpsNeg","Epithelial Cells_LpsPos","Fras1+ Epithelial Cells_LpsNeg","Fras1+ Epithelial Cells_LpsPos",
                                                                                        "Endothelial Cells_LpsNeg","Endothelial Cells_LpsPos","Vascular Associated Cells_LpsNeg","Vascular Associated Cells_LpsPos",
                                                                                        "Stromal Fibroblasts_LpsNeg","Stromal Fibroblasts_LpsPos","Stalk Fibroblasts_LpsNeg","Stalk Fibroblasts_LpsPos",
                                                                                        "Fibroblast-like Cells_LpsNeg","Fibroblast-like Cells_LpsPos","Macrophages_LpsNeg","Macrophages_LpsPos","Epiplexus Cells_LpsNeg", 
                                                                                        "NK Cells_LpsNeg","NK Cells_LpsPos","Dendritic Cells_LpsNeg","Dendritic Cells_LpsPos",
                                                                                        "T Cells_LpsNeg","T Cells_LpsPos","Monocytes_LpsNeg","Monocytes_LpsPos","Neutrophils_LpsPos","Mitotic Cells_LpsNeg","Mitotic Cells_LpsPos",
                                                                                        "Neuronal and Glial Cells_LpsNeg","Neuronal and Glial Cells_LpsPos"))
## Now with numbers
seuratObjNew2<-seuratObjNew
seuratObjNew2@meta.data$newAnnotated_clusters<-seuratObjNew2@meta.data$newClusters #Save split annotation separately!
levels(seuratObjNew2@meta.data$newClusters)<-seq(1,30) #Numbered annotation!!
Idents(seuratObjNew2)<-seuratObjNew2@meta.data$newClusters

## Reorganize annotated clusters too
levels(seuratObjNew2$annotated_clusters)
seuratObjNew2$annotated_clusters<-factor(as.character(seuratObjNew2$annotated_clusters), levels = c("Epithelial Cells","Fras1+ Epithelial Cells","Endothelial Cells","Vascular Associated Cells",
                                                                                                    "Stromal Fibroblasts","Stalk Fibroblasts","Fibroblast-like Cells","Macrophages","Epiplexus Cells", 
                                                                                                    "NK Cells","Dendritic Cells","T Cells","Monocytes","Neutrophils","Mitotic Cells","Neuronal and Glial Cells"))

## Visualize
Colorset<-c(brewer.pal(12,"Paired"),"#d3d3d3","#808080","magenta","magenta4","Black","#00ffff","turquoise4","gold","goldenrod2",
            "Olivedrab1", "Olivedrab3", "Orangered", "Darkred","Blue4","Hotpink1","Hotpink3", "Seagreen1", "Seagreen4")

pdf(file=paste0(sampleFolder,"results/QC/Paper/15_Numbered_UMAP_paper_110620.pdf"), width = 17, height = 12)
DimPlot(seuratObjNew2, reduction = "umap", label = T, repel = T, label.size = 8, cols = Colorset)
dev.off()

## Frequency bar graph for paper to make pdf
seuratObjNew2@meta.data$Treatment<-factor(seuratObjNew2@meta.data$Treatment, levels=c("LpsPos","LpsNeg"))
Condition <- seuratObjNew2@meta.data$Treatment
cluster <- seuratObjNew2$annotated_clusters

data <- data.frame(table(Condition, cluster))

col_freq<-c("#EE5252","#C0BEBE") #030720

pdf(file=paste0(sampleFolder,"results/QC/Paper/14_SampleDistribution_ggplot2_paper_new_flipped.pdf"), width = 8, height = 8)
ggplot(data, aes(fill=Condition, y=Freq, x=cluster)) + theme_bw() + coord_flip() + scale_x_discrete(limits=rev) +
  #scale_x_continuous(position = "top") +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(angle = 90, hjust = 0.5)) +
  geom_bar(position="fill", stat="identity", colour="white") +
  scale_fill_manual(values=col_freq)
dev.off()

# Revert back factor
seuratObjNew2@meta.data$Treatment<-factor(seuratObjNew2@meta.data$Treatment, levels=c("LpsNeg","LpsPos"))

######################3

### Dotplot paper
# CPE: "Otx2"
# Fras1+ CPE: "Fras1"
# EC: "Plvap"
# VAC: "Kcnj8"
# FB I: "Dpep1" (stromal)
# FB II: "Igfbp6" (stalk)
# FB like: "Smoc2"
# MF: "Adgre1" 
# Epiplexus cells: "P2ry12"
# NK: "Klrb1c"
# DCs: "Cd209a","Ccr7","Xcr1"
# T_cells: "Icos"
# MC: "Ly6c2"
# NF: "S100a8"
# Mitotic: "Mki67"
# Neuronal and glial: "Tubb3" + "Fabp7"

## Dotplots
Colors_dotplot<-c("#071AE5","#F50635") 

wantedGenes<-c("Otx2","Fras1","Plvap","Kcnj8","Dpep1","Igfbp6","Smoc2","Adgre1","P2ry12",
               "Klrb1c","Cd209a","Ccr7","Xcr1","Icos","Ly6c2","S100a8","Mki67","Tubb3","Fabp7")
wantedGenes<-rev(wantedGenes)

D1<-DotPlot(seuratObjNew2, features = wantedGenes, cols = Colors_dotplot, dot.scale = 12)
D2<-ggpar(D1, orientation = "horizontal") #ggpubr package to reorient ggplot objects

pdf(file=paste0(sampleFolder,"results/QC/Paper/16_Dotplot_Otx2_horizontal_030720.pdf"), width = 15, height = 15)
D2
dev.off()

#####################################################

## Final annotation update (07/2023)
# Fibroblast-like Cells → ABCs
# Stalk Fibroblasts → BBCs
# Fras1+ Epithelial cells → ??? (not sure what else to call. Should not be removed and should remain separate from CPE. Leave it as is!)
# Microglia-like Macrophages → Epiplexus Cells
levels(seuratObjNew2$annotated_clusters)[6]<-"Base Barrier Cells"
levels(seuratObjNew2$annotated_clusters)[7]<-"Arachnoid Barrier Cells"
levels(seuratObjNew2$newAnnotated_clusters)[11]<-"Base Barrier Cells_LpsNeg"
levels(seuratObjNew2$newAnnotated_clusters)[12]<-"Base Barrier Cells_LpsPos"
levels(seuratObjNew2$newAnnotated_clusters)[13]<-"Arachnoid Barrier Cells_LpsNeg"
levels(seuratObjNew2$newAnnotated_clusters)[14]<-"Arachnoid Barrier Cells_LpsPos"

########################################
##### Paper marker analysis
########################################

## New marker analysis with final annotation (split and not split per condition!)
### Find markers for every cluster compared to all remaining cells

## First round
Idents(seuratObjNew2)<-seuratObjNew2$annotated_clusters
allMarkers_v1 <- FindAllMarkers(seuratObjNew2, assay = "RNA", only.pos = T)
table(allMarkers_v1$cluster)
saveRDS(allMarkers_v1, file=paste0(sampleFolder,"Robjects/allMarkers_annotated_v1_",sampleName,".rds"))

### Create list with markers
markersList_v1<-list()

for(i in names(table(allMarkers_v1$cluster))){
  tmp<-allMarkers_v1[allMarkers_v1$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  markersList_v1[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

### Write to Excel
write.xlsx(markersList_v1, file =paste0(sampleFolder, "results/allClusters_annotated_v1_",sampleName,".xlsx"))

## Second round
Idents(seuratObjNew2)<-seuratObjNew2$newAnnotated_clusters
allMarkers_v2 <- FindAllMarkers(seuratObjNew2, assay = "RNA", only.pos = T)
table(allMarkers_v2$cluster)
saveRDS(allMarkers_v2, file=paste0(sampleFolder,"Robjects/allMarkers_annotated_v2_split_",sampleName,".rds"))

### Create list with markers
markersList_v2<-list()

for(i in names(table(allMarkers_v2$cluster))){
  tmp<-allMarkers_v2[allMarkers_v2$cluster==i,]
  tmp$score<-tmp$pct.1/(tmp$pct.2+0.01)*tmp$avg_logFC
  
  markersList_v2[[i]]<-tmp[order(tmp$score, decreasing=TRUE),]
}

## Adjust names (too long for excel tabs)
for (i in names(markersList_v2)){
  print(nchar(i))
} #7/8/

names(markersList_v2)[7]<-"VAC_LpsNeg"
names(markersList_v2)[8]<-"VAC_LpsPos"

### Write to Excel
write.xlsx(markersList_v2, file =paste0(sampleFolder, "results/allClusters_annotated_v2_split_",sampleName,".xlsx"))

#####################################################

##### Save object
saveRDS(seuratObjNew2, file=paste0(sampleFolder,"Robjects/seuratObj_subset_paper_numbered_",sampleName,".rds"))

##### Read object
seuratObjNew2<-readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_subset_paper_numbered_",sampleName,".rds"))

################################################################################################################
################################################################################################################

## Create diet object for online tool 
seuratObj_subset_diet<-DietSeurat(seuratObjNew2, counts = T, data = T, scale.data = F,
                                  assays = c("RNA"), dimreducs = "umap", graphs = NULL)

## New metadata names
seuratObj_subset_diet$Sample<-as.factor(seuratObj_subset_diet$orig.ident)
seuratObj_subset_diet$Treatment
seuratObj_subset_diet$Annotation<-seuratObj_subset_diet$annotated_clusters
seuratObj_subset_diet$Annotation_split<-seuratObj_subset_diet$newAnnotated_clusters
seuratObj_subset_diet$Annotation_numbered<-seuratObj_subset_diet$newClusters

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

col2hex("Seagreen4")
All<-c(levels(seuratObj_subset_diet$Annotation_split),levels(seuratObj_subset_diet$Annotation_numbered),
       levels(seuratObj_subset_diet$Annotation),levels(seuratObj_subset_diet$Treatment),
       levels(seuratObj_subset_diet$Sample))
Color_info<-c(rep(c(brewer.pal(12,"Paired"),"#d3d3d3","#808080","#FF00FF","#8B008B","#000000","#00ffff","#00868B","#FFD700","#EEB422",
                "#C0FF3E", "#9ACD32", "#FF4500", "#8B0000","#00008B","#FF6EB4","#CD6090", "#54FF9F", "#2E8B57"),2),
              gg_color_hue(16), c("#C0BEBE","#EE5252"), gg_color_hue(4))
Metadata_column<-c(rep("Annotation_split",30),rep("Annotation_numbered",30),
                   rep("Annotation",16),rep("Treatment",2),
                   rep("Sample",4))
Info_Kevin<-as.data.frame(cbind(All,Color_info,Metadata_column))

write.xlsx(Info_Kevin, file =paste0(sampleFolder, "results/QC/Info_Kevin_object_",sampleName,".xlsx"))

##### Save object
saveRDS(seuratObj_subset_diet, file=paste0(sampleFolder,"Robjects/seuratObj_subset_paper_diet_",sampleName,"_2023.rds"))

##### Read object
seuratObj_subset_diet<-readRDS(file=paste0(sampleFolder,"Robjects/seuratObj_paper_diet_",sampleName,"_2023.rds"))

##########################################################################################

########################################
##### Differential markers
########################################

### Still using the original annotation (not the updated names!)
Idents(seuratObjNew2)<-seuratObjNew2@meta.data$newClusters

########## 2. GET MARKERS ########## 
getDEgenes<-function(ident1, ident2){
  markersDiff <- FindMarkers(seuratObjNew2, ident.1 = ident1, ident.2 = ident2,
                             min.pct = 0.10) # No min diff pct! logFC default threshold is 0.25
  markersDiff<-markersDiff[markersDiff$p_val_adj < 0.01,]
  markersDiff<-markersDiff[order(markersDiff$avg_logFC, decreasing = T),]

  markersDiff$geneSymbol<-rownames(markersDiff)
  markersDiff$pct.1<-markersDiff$pct.1+0.001
  markersDiff$pct.2<-markersDiff$pct.2+0.001

  markersDiff<-rbind(markersDiff[markersDiff$avg_logFC > 0,] %>% dplyr::mutate(.,score=pct.1/pct.2*avg_logFC),
                     markersDiff[markersDiff$avg_logFC < 0,] %>% dplyr::mutate(.,score=pct.2/pct.1*avg_logFC))
  markersDiff<-markersDiff[order(markersDiff$score, decreasing = T),]
  return(markersDiff)
}

#### Get diff markers LpsPos vs LpsNeg #####
Epithelialcells_LpsPosvsLpsNeg<-getDEgenes("Epithelial Cells_LpsPos","Epithelial Cells_LpsNeg")
Epithelialcells_LpsPosvsLpsNeg<-Epithelialcells_LpsPosvsLpsNeg[order(Epithelialcells_LpsPosvsLpsNeg$avg_logFC,decreasing = T),]
head(Epithelialcells_LpsPosvsLpsNeg)
dim(Epithelialcells_LpsPosvsLpsNeg)

Macrophages_LpsPosvsLpsNeg<-getDEgenes("Macrophages_LpsPos","Macrophages_LpsNeg")
Macrophages_LpsPosvsLpsNeg<-Macrophages_LpsPosvsLpsNeg[order(Macrophages_LpsPosvsLpsNeg$avg_logFC,decreasing = T),]
head(Macrophages_LpsPosvsLpsNeg)
dim(Macrophages_LpsPosvsLpsNeg)

Endothelialcells_LpsPosvsLpsNeg<-getDEgenes("Endothelial Cells_LpsPos","Endothelial Cells_LpsNeg")
Endothelialcells_LpsPosvsLpsNeg<-Endothelialcells_LpsPosvsLpsNeg[order(Endothelialcells_LpsPosvsLpsNeg$avg_logFC,decreasing = T),]
head(Endothelialcells_LpsPosvsLpsNeg)
dim(Endothelialcells_LpsPosvsLpsNeg)

Fibroblasts_Type_1_LpsPosvsLpsNeg<-getDEgenes("Fibroblasts Type 1_LpsPos","Fibroblasts Type 1_LpsNeg")
Fibroblasts_Type_1_LpsPosvsLpsNeg<-Fibroblasts_Type_1_LpsPosvsLpsNeg[order(Fibroblasts_Type_1_LpsPosvsLpsNeg$avg_logFC,decreasing = T),]
head(Fibroblasts_Type_1_LpsPosvsLpsNeg)
dim(Fibroblasts_Type_1_LpsPosvsLpsNeg)

Fibroblasts_Type_2_LpsPosvsLpsNeg<-getDEgenes("Fibroblasts Type 2_LpsPos","Fibroblasts Type 2_LpsNeg")
Fibroblasts_Type_2_LpsPosvsLpsNeg<-Fibroblasts_Type_2_LpsPosvsLpsNeg[order(Fibroblasts_Type_2_LpsPosvsLpsNeg$avg_logFC,decreasing = T),]
head(Fibroblasts_Type_2_LpsPosvsLpsNeg)
dim(Fibroblasts_Type_2_LpsPosvsLpsNeg)

Vascular_associated_cells_LpsPosvsLpsNeg<-getDEgenes("Vascular Associated Cells_LpsPos","Vascular Associated Cells_LpsNeg")
Vascular_associated_cells_LpsPosvsLpsNeg<-Vascular_associated_cells_LpsPosvsLpsNeg[order(Vascular_associated_cells_LpsPosvsLpsNeg$avg_logFC,decreasing = T),]
head(Vascular_associated_cells_LpsPosvsLpsNeg)
dim(Vascular_associated_cells_LpsPosvsLpsNeg)

NK_cells_LpsPosvsLpsNeg<-getDEgenes("NK Cells_LpsPos","NK Cells_LpsNeg")
NK_cells_LpsPosvsLpsNeg<-NK_cells_LpsPosvsLpsNeg[order(NK_cells_LpsPosvsLpsNeg$avg_logFC,decreasing = T),]
head(NK_cells_LpsPosvsLpsNeg)
dim(NK_cells_LpsPosvsLpsNeg)

Dendritic_Cells_LpsPosvsLpsNeg<-getDEgenes("Dendritic Cells_LpsPos","Dendritic Cells_LpsNeg")
Dendritic_Cells_LpsPosvsLpsNeg<-Dendritic_Cells_LpsPosvsLpsNeg[order(Dendritic_Cells_LpsPosvsLpsNeg$avg_logFC,decreasing = T),]
head(Dendritic_Cells_LpsPosvsLpsNeg)
dim(Dendritic_Cells_LpsPosvsLpsNeg)

Mitotic_Cells_LpsPosvsLpsNeg<-getDEgenes("Mitotic Cells_LpsPos","Mitotic Cells_LpsNeg")
Mitotic_Cells_LpsPosvsLpsNeg<-Mitotic_Cells_LpsPosvsLpsNeg[order(Mitotic_Cells_LpsPosvsLpsNeg$avg_logFC,decreasing = T),]
head(Mitotic_Cells_LpsPosvsLpsNeg)
dim(Mitotic_Cells_LpsPosvsLpsNeg)

T_Cells_LpsPosvsLpsNeg<-getDEgenes("T Cells_LpsPos","T Cells_LpsNeg")
T_Cells_LpsPosvsLpsNeg<-T_Cells_LpsPosvsLpsNeg[order(T_Cells_LpsPosvsLpsNeg$avg_logFC,decreasing = T),]
head(T_Cells_LpsPosvsLpsNeg)
dim(T_Cells_LpsPosvsLpsNeg)

##add to list
listDiffMarkers<-tibble::lst(Epithelialcells_LpsPosvsLpsNeg, Macrophages_LpsPosvsLpsNeg,Endothelialcells_LpsPosvsLpsNeg, 
                             Fibroblasts_Type_1_LpsPosvsLpsNeg, Fibroblasts_Type_2_LpsPosvsLpsNeg, Vascular_associated_cells_LpsPosvsLpsNeg,
                             NK_cells_LpsPosvsLpsNeg, Dendritic_Cells_LpsPosvsLpsNeg,
                             Mitotic_Cells_LpsPosvsLpsNeg,T_Cells_LpsPosvsLpsNeg)
lapply(listDiffMarkers, dim)
listDiffMarkers<-lapply(listDiffMarkers,function(x){x<-x[order(x$score, decreasing=T),]})

#Check settings
saveRDS(listDiffMarkers,file=paste0(sampleFolder,"Robjects/markersDiffSamples_Full.rds"))

### Write to Excel
write.xlsx(listDiffMarkers, file = paste0(sampleFolder,"results/summaryDiffMarkers_Full.xlsx"))

###################

#### Repeat but now split per ventricle (also with original annotation, not the updated names!)
seuratObjNew2@meta.data$newClusters2<-paste0(seuratObjNew2@meta.data$annotated_clusters,"_",seuratObjNew2@meta.data$orig.ident)
Idents(seuratObjNew2)<-seuratObjNew2@meta.data$newClusters2

#### Get diff markers LpsPosLV vs LpsNegLV #####
Epithelialcells_LpsPosLVvsLpsNegLV<-getDEgenes("Epithelial Cells_LpsPosLat","Epithelial Cells_LpsNegLat")
Epithelialcells_LpsPosLVvsLpsNegLV<-Epithelialcells_LpsPosLVvsLpsNegLV[order(Epithelialcells_LpsPosLVvsLpsNegLV$avg_logFC,decreasing = T),]
head(Epithelialcells_LpsPosLVvsLpsNegLV)
dim(Epithelialcells_LpsPosLVvsLpsNegLV)

Xist_CPE_LpsPosLVvsLpsNegLV<-getDEgenes("Xist+ Epithelial Cells_LpsPosLat","Xist+ Epithelial Cells_LpsNegLat")
Xist_CPE_LpsPosLVvsLpsNegLV<-Xist_CPE_LpsPosLVvsLpsNegLV[order(Xist_CPE_LpsPosLVvsLpsNegLV$avg_logFC,decreasing = T),]
head(Xist_CPE_LpsPosLVvsLpsNegLV)
dim(Xist_CPE_LpsPosLVvsLpsNegLV)

Macrophages_LpsPosLVvsLpsNegLV<-getDEgenes("Macrophages_LpsPosLat","Macrophages_LpsNegLat")
Macrophages_LpsPosLVvsLpsNegLV<-Macrophages_LpsPosLVvsLpsNegLV[order(Macrophages_LpsPosLVvsLpsNegLV$avg_logFC,decreasing = T),]
head(Macrophages_LpsPosLVvsLpsNegLV)
dim(Macrophages_LpsPosLVvsLpsNegLV)

Endothelialcells_LpsPosLVvsLpsNegLV<-getDEgenes("Endothelial Cells_LpsPosLat","Endothelial Cells_LpsNegLat")
Endothelialcells_LpsPosLVvsLpsNegLV<-Endothelialcells_LpsPosLVvsLpsNegLV[order(Endothelialcells_LpsPosLVvsLpsNegLV$avg_logFC,decreasing = T),]
head(Endothelialcells_LpsPosLVvsLpsNegLV)
dim(Endothelialcells_LpsPosLVvsLpsNegLV)

Fibroblasts_Type_1_LpsPosLVvsLpsNegLV<-getDEgenes("Fibroblasts Type 1_LpsPosLat","Fibroblasts Type 1_LpsNegLat")
Fibroblasts_Type_1_LpsPosLVvsLpsNegLV<-Fibroblasts_Type_1_LpsPosLVvsLpsNegLV[order(Fibroblasts_Type_1_LpsPosLVvsLpsNegLV$avg_logFC,decreasing = T),]
head(Fibroblasts_Type_1_LpsPosLVvsLpsNegLV)
dim(Fibroblasts_Type_1_LpsPosLVvsLpsNegLV)

Vascular_associated_cells_LpsPosLVvsLpsNegLV<-getDEgenes("Vascular Associated Cells_LpsPosLat","Vascular Associated Cells_LpsNegLat")
Vascular_associated_cells_LpsPosLVvsLpsNegLV<-Vascular_associated_cells_LpsPosLVvsLpsNegLV[order(Vascular_associated_cells_LpsPosLVvsLpsNegLV$avg_logFC,decreasing = T),]
head(Vascular_associated_cells_LpsPosLVvsLpsNegLV)
dim(Vascular_associated_cells_LpsPosLVvsLpsNegLV)

NK_cells_LpsPosLVvsLpsNegLV<-getDEgenes("NK Cells_LpsPosLat","NK Cells_LpsNegLat")
NK_cells_LpsPosLVvsLpsNegLV<-NK_cells_LpsPosLVvsLpsNegLV[order(NK_cells_LpsPosLVvsLpsNegLV$avg_logFC,decreasing = T),]
head(NK_cells_LpsPosLVvsLpsNegLV)
dim(NK_cells_LpsPosLVvsLpsNegLV)

Dendritic_Cells_LpsPosLVvsLpsNegLV<-getDEgenes("Dendritic Cells_LpsPosLat","Dendritic Cells_LpsNegLat")
Dendritic_Cells_LpsPosLVvsLpsNegLV<-Dendritic_Cells_LpsPosLVvsLpsNegLV[order(Dendritic_Cells_LpsPosLVvsLpsNegLV$avg_logFC,decreasing = T),]
head(Dendritic_Cells_LpsPosLVvsLpsNegLV)
dim(Dendritic_Cells_LpsPosLVvsLpsNegLV)

Mitotic_Cells_LpsPosLVvsLpsNegLV<-getDEgenes("Mitotic Cells_LpsPosLat","Mitotic Cells_LpsNegLat")
Mitotic_Cells_LpsPosLVvsLpsNegLV<-Mitotic_Cells_LpsPosLVvsLpsNegLV[order(Mitotic_Cells_LpsPosLVvsLpsNegLV$avg_logFC,decreasing = T),]
head(Mitotic_Cells_LpsPosLVvsLpsNegLV)
dim(Mitotic_Cells_LpsPosLVvsLpsNegLV)

T_Cells_LpsPosLVvsLpsNegLV<-getDEgenes("T Cells_LpsPosLat","T Cells_LpsNegLat")
T_Cells_LpsPosLVvsLpsNegLV<-T_Cells_LpsPosLVvsLpsNegLV[order(T_Cells_LpsPosLVvsLpsNegLV$avg_logFC,decreasing = T),]
head(T_Cells_LpsPosLVvsLpsNegLV)
dim(T_Cells_LpsPosLVvsLpsNegLV)

#### Get diff markers LpsPos4V vs LpsNeg4V #####
Epithelialcells_LpsPos4VvsLpsNeg4V<-getDEgenes("Epithelial Cells_LpsPosFour","Epithelial Cells_LpsNegFour")
Epithelialcells_LpsPos4VvsLpsNeg4V<-Epithelialcells_LpsPos4VvsLpsNeg4V[order(Epithelialcells_LpsPos4VvsLpsNeg4V$avg_logFC,decreasing = T),]
head(Epithelialcells_LpsPos4VvsLpsNeg4V)
dim(Epithelialcells_LpsPos4VvsLpsNeg4V)

Xist_CPE_LpsPos4VvsLpsNeg4V<-getDEgenes("Xist+ Epithelial Cells_LpsPosFour","Xist+ Epithelial Cells_LpsNegFour")
Xist_CPE_LpsPos4VvsLpsNeg4V<-Xist_CPE_LpsPos4VvsLpsNeg4V[order(Xist_CPE_LpsPos4VvsLpsNeg4V$avg_logFC,decreasing = T),]
head(Xist_CPE_LpsPos4VvsLpsNeg4V)
dim(Xist_CPE_LpsPos4VvsLpsNeg4V)

Macrophages_LpsPos4VvsLpsNeg4V<-getDEgenes("Macrophages_LpsPosFour","Macrophages_LpsNegFour")
Macrophages_LpsPos4VvsLpsNeg4V<-Macrophages_LpsPos4VvsLpsNeg4V[order(Macrophages_LpsPos4VvsLpsNeg4V$avg_logFC,decreasing = T),]
head(Macrophages_LpsPos4VvsLpsNeg4V)
dim(Macrophages_LpsPos4VvsLpsNeg4V)

Endothelialcells_LpsPos4VvsLpsNeg4V<-getDEgenes("Endothelial Cells_LpsPosFour","Endothelial Cells_LpsNegFour")
Endothelialcells_LpsPos4VvsLpsNeg4V<-Endothelialcells_LpsPos4VvsLpsNeg4V[order(Endothelialcells_LpsPos4VvsLpsNeg4V$avg_logFC,decreasing = T),]
head(Endothelialcells_LpsPos4VvsLpsNeg4V)
dim(Endothelialcells_LpsPos4VvsLpsNeg4V)

Fibroblasts_Type_1_LpsPos4VvsLpsNeg4V<-getDEgenes("Fibroblasts Type 1_LpsPosFour","Fibroblasts Type 1_LpsNegFour")
Fibroblasts_Type_1_LpsPos4VvsLpsNeg4V<-Fibroblasts_Type_1_LpsPos4VvsLpsNeg4V[order(Fibroblasts_Type_1_LpsPos4VvsLpsNeg4V$avg_logFC,decreasing = T),]
head(Fibroblasts_Type_1_LpsPos4VvsLpsNeg4V)
dim(Fibroblasts_Type_1_LpsPos4VvsLpsNeg4V)

Fibroblasts_Type_2_LpsPos4VvsLpsNeg4V<-getDEgenes("Fibroblasts Type 2_LpsPosFour","Fibroblasts Type 2_LpsNegFour")
Fibroblasts_Type_2_LpsPos4VvsLpsNeg4V<-Fibroblasts_Type_2_LpsPos4VvsLpsNeg4V[order(Fibroblasts_Type_2_LpsPos4VvsLpsNeg4V$avg_logFC,decreasing = T),]
head(Fibroblasts_Type_2_LpsPos4VvsLpsNeg4V)
dim(Fibroblasts_Type_2_LpsPos4VvsLpsNeg4V)

Vascular_associated_cells_LpsPos4VvsLpsNeg4V<-getDEgenes("Vascular Associated Cells_LpsPosFour","Vascular Associated Cells_LpsNegFour")
Vascular_associated_cells_LpsPos4VvsLpsNeg4V<-Vascular_associated_cells_LpsPos4VvsLpsNeg4V[order(Vascular_associated_cells_LpsPos4VvsLpsNeg4V$avg_logFC,decreasing = T),]
head(Vascular_associated_cells_LpsPos4VvsLpsNeg4V)
dim(Vascular_associated_cells_LpsPos4VvsLpsNeg4V)

NK_cells_LpsPos4VvsLpsNeg4V<-getDEgenes("NK Cells_LpsPosFour","NK Cells_LpsNegFour")
NK_cells_LpsPos4VvsLpsNeg4V<-NK_cells_LpsPos4VvsLpsNeg4V[order(NK_cells_LpsPos4VvsLpsNeg4V$avg_logFC,decreasing = T),]
head(NK_cells_LpsPos4VvsLpsNeg4V)
dim(NK_cells_LpsPos4VvsLpsNeg4V)

Dendritic_Cells_LpsPos4VvsLpsNeg4V<-getDEgenes("Dendritic Cells_LpsPosFour","Dendritic Cells_LpsNegFour")
Dendritic_Cells_LpsPos4VvsLpsNeg4V<-Dendritic_Cells_LpsPos4VvsLpsNeg4V[order(Dendritic_Cells_LpsPos4VvsLpsNeg4V$avg_logFC,decreasing = T),]
head(Dendritic_Cells_LpsPos4VvsLpsNeg4V)
dim(Dendritic_Cells_LpsPos4VvsLpsNeg4V)

Mitotic_Cells_LpsPos4VvsLpsNeg4V<-getDEgenes("Mitotic Cells_LpsPosFour","Mitotic Cells_LpsNegFour")
Mitotic_Cells_LpsPos4VvsLpsNeg4V<-Mitotic_Cells_LpsPos4VvsLpsNeg4V[order(Mitotic_Cells_LpsPos4VvsLpsNeg4V$avg_logFC,decreasing = T),]
head(Mitotic_Cells_LpsPos4VvsLpsNeg4V)
dim(Mitotic_Cells_LpsPos4VvsLpsNeg4V)

T_Cells_LpsPos4VvsLpsNeg4V<-getDEgenes("T Cells_LpsPosFour","T Cells_LpsNegFour")
T_Cells_LpsPos4VvsLpsNeg4V<-T_Cells_LpsPos4VvsLpsNeg4V[order(T_Cells_LpsPos4VvsLpsNeg4V$avg_logFC,decreasing = T),]
head(T_Cells_LpsPos4VvsLpsNeg4V)
dim(T_Cells_LpsPos4VvsLpsNeg4V)


##add to list
listDiffMarkers_4V_LV<-tibble::lst(Epithelialcells_LpsPosLVvsLpsNegLV, Xist_CPE_LpsPosLVvsLpsNegLV,
                                   Macrophages_LpsPosLVvsLpsNegLV,Endothelialcells_LpsPosLVvsLpsNegLV, 
                                   Fibroblasts_Type_1_LpsPosLVvsLpsNegLV, Vascular_associated_cells_LpsPosLVvsLpsNegLV,
                                   NK_cells_LpsPosLVvsLpsNegLV, Dendritic_Cells_LpsPosLVvsLpsNegLV,
                                   Epithelialcells_LpsPos4VvsLpsNeg4V, Xist_CPE_LpsPos4VvsLpsNeg4V,
                                   Macrophages_LpsPos4VvsLpsNeg4V,Endothelialcells_LpsPos4VvsLpsNeg4V, 
                                   Fibroblasts_Type_1_LpsPos4VvsLpsNeg4V, Fibroblasts_Type_2_LpsPos4VvsLpsNeg4V, Vascular_associated_cells_LpsPos4VvsLpsNeg4V,
                                   NK_cells_LpsPos4VvsLpsNeg4V, Dendritic_Cells_LpsPos4VvsLpsNeg4V)

lapply(listDiffMarkers_4V_LV, dim)
listDiffMarkers_4V_LV<-lapply(listDiffMarkers_4V_LV,function(x){x<-x[order(x$score, decreasing=T),]})

saveRDS(listDiffMarkers_4V_LV,file=paste0(sampleFolder,"Robjects/markersDiffSamples_4V_LV_Full.rds"))

### Write to Excel
library('openxlsx')
write.xlsx(listDiffMarkers_4V_LV, file = paste0(sampleFolder,"results/summaryDiffMarkers_4V_LV_Full.xlsx"))

###################################################

## Combine the differential markers into 1 tab per cell type

listDiffMarkers_combined<-tibble::lst()

Cell_pops<-c('Epithelialcells', "Xist_CPE", "Macrophages", "Endothelialcells", "Fibroblasts_Type_1",
             "Fibroblasts_Type_2", "Vascular_associated_cells", "NK_cells", "Dendritic_Cells")

for (pop in Cell_pops){
  indices<- grep(pop, names(listDiffMarkers_4V_LV)) #Check for cell populations (in both lists or not?)
  print(indices)
  if (length(indices)==2){ #In both lists
    listDiffMarkers_4V_LV[[indices[1]]]$Cell_pop<-pop
    listDiffMarkers_4V_LV[[indices[2]]]$Cell_pop<-pop
    listDiffMarkers_combined[[pop]]<-rbind(listDiffMarkers_4V_LV[[indices[1]]],listDiffMarkers_4V_LV[[indices[2]]]) #Combine the lists into 1
  } else { #Only one list
    listDiffMarkers_4V_LV[[indices[1]]]$Cell_pop<-pop
    listDiffMarkers_combined[[pop]]<-listDiffMarkers_4V_LV[[indices[1]]]
  }
}

for (y in 1:length(listDiffMarkers_combined)) {
  genes_pop<-unique(listDiffMarkers_combined[[y]]$geneSymbol) #Unique list of genes in cell population
  listDiffMarkers_combined[[y]]$Ventricle_unique<-F #Start off with all False
  for (gene in genes_pop) {
    gene<-paste0("\\<",gene,"\\>") #Direct match needed!!
    indices_gene<-grep(gene,listDiffMarkers_combined[[y]]$geneSymbol)
    if (length(indices_gene)==1) { #If only one match found, then unique
      listDiffMarkers_combined[[y]]$Ventricle_unique[indices_gene[1]]<-T #Update column
    }
  }
}

#Combine all lists into one
listDiffMarkers_all<-rbind(listDiffMarkers_combined[[1]],listDiffMarkers_combined[[2]],listDiffMarkers_combined[[3]],
                           listDiffMarkers_combined[[4]],listDiffMarkers_combined[[5]],listDiffMarkers_combined[[6]],
                           listDiffMarkers_combined[[7]],listDiffMarkers_combined[[8]],listDiffMarkers_combined[[9]])

# Get the full list of unique genes (over all cell populations)
gene_list<-listDiffMarkers_all$geneSymbol
gene_list_unique<-unique(gene_list)
gene_df<-data.frame()

#Perform check for strange symbols in genes
gene_list_unique[grepl('\\(', gene_list_unique)]

for (gene in 1:length(gene_list_unique)){ #Count in how many lists the gene appears
  count_gene<-0 #Reset for each gene
  
  for (y in 1:length(listDiffMarkers_combined)) { #Check in each list
    genes_pop<-unique(listDiffMarkers_combined[[y]]$geneSymbol) #Get unique list in population
    if(!grepl('\\(', gene_list_unique[gene])){ ##To avoid issue with one gene (ROSA)!!
      gene_to_search<-paste0("\\<",gene_list_unique[gene],"\\>")
      if(length(grep(gene_to_search,genes_pop))==1){
        count_gene<-count_gene+1  #Add to count
      }
    } else{ #Gene which gives issue
      gene_to_search<-"ROSA"
      if(length(grep(gene_to_search,genes_pop))==1){
        count_gene<-count_gene+1
      }
    }

  }
  
  gene_df[gene,"Gene_name"]<-gene_list_unique[gene] #Add name to dataframe
  gene_df[gene,"Gene_count"]<-count_gene #Add final count to dataframe
}

for (y in 1:length(listDiffMarkers_combined)) { #Go over each list
  listDiffMarkers_combined[[y]]$Gene_pop_count<-1 #Default value
  listDiffMarkers_combined[[y]]$Shared_gene<-F #Default value
  listDiffMarkers_combined[[y]]$Cell_pop_Unique<-F #Default value
  listDiffMarkers_combined[[y]]$Gene_class<-1 #Default value
  #4 types of genes: 
  #Basic genes: class 1
  #Shared gene: testData[,5]==T -> class 2
  #Not shared, cell type specific: testData[,6]==T -> class 3
  #Cell type and ventricle specific: testData[,4]==T & testData[,6]==T -> class 4
  
  
  for (x in 1:nrow(listDiffMarkers_combined[[y]])) { #Go over each row and update values
    if(!grepl('\\(', listDiffMarkers_combined[[y]]$geneSymbol[x])){ ##To avoid issue with one gene!!
      gene<-paste0("\\<",listDiffMarkers_combined[[y]]$geneSymbol[x],"\\>") #Need perfect match
      if(grep(gene,gene_df[,1])){ #Look for gene in gene dataframe
        listDiffMarkers_combined[[y]]$Gene_pop_count[x]<-gene_df[grep(gene,gene_df[,1]),2] #Update Gene_pop_count value with value dataframe
        
        if(gene_df[grep(gene,gene_df[,1]),2]==1){ #If that value is 1, then update column Cell_pop_unique
          listDiffMarkers_combined[[y]]$Cell_pop_Unique[x]<-T
        } else if (gene_df[grep(gene,gene_df[,1]),2]>3){ #If that value is larger than 3, then update column Shared_gene
          listDiffMarkers_combined[[y]]$Shared_gene[x]<-T
        }
        
        if(listDiffMarkers_combined[[y]]$Shared_gene[x]){ #If Shared_gene Column is true, then update Gene_class value to 2
          listDiffMarkers_combined[[y]]$Gene_class[x]<-2
        } else if (listDiffMarkers_combined[[y]]$Cell_pop_Unique[x]){ #If Cell_pop_unique Column is true, then update Gene_class to 3
          listDiffMarkers_combined[[y]]$Gene_class[x]<-3
        } 
        #separate if necessary!!
        if ((listDiffMarkers_combined[[y]]$Ventricle_unique[x]) && (listDiffMarkers_combined[[y]]$Cell_pop_Unique[x])){ #Also vent. unique!
          listDiffMarkers_combined[[y]]$Gene_class[x]<-4 #Update Gene_class to 4
        }
      }
    } else { # Repeat for issue gene!
      gene<-"ROSA"
      if(grep(gene,gene_df[,1])){
        listDiffMarkers_combined[[y]]$Gene_pop_count[x]<-gene_df[grep(gene,gene_df[,1]),2]
        
        if(gene_df[grep(gene,gene_df[,1]),2]==1){
          listDiffMarkers_combined[[y]]$Cell_pop_Unique[x]<-T
        } else if (gene_df[grep(gene,gene_df[,1]),2]>3){
          listDiffMarkers_combined[[y]]$Shared_gene[x]<-T
        }
        
        if(listDiffMarkers_combined[[y]]$Shared_gene[x]){
          listDiffMarkers_combined[[y]]$Gene_class[x]<-2
        } else if (listDiffMarkers_combined[[y]]$Cell_pop_Unique[x]){
          listDiffMarkers_combined[[y]]$Gene_class[x]<-3
        } 
        
        if ((listDiffMarkers_combined[[y]]$Ventricle_unique[x]) && (listDiffMarkers_combined[[y]]$Cell_pop_Unique[x])){
          listDiffMarkers_combined[[y]]$Gene_class[x]<-4
        }
      }
    }
  }
}

# Remake full list with new info added
listDiffMarkers_all<-rbind(listDiffMarkers_combined[[1]],listDiffMarkers_combined[[2]],listDiffMarkers_combined[[3]],
                           listDiffMarkers_combined[[4]],listDiffMarkers_combined[[5]],listDiffMarkers_combined[[6]],
                           listDiffMarkers_combined[[7]],listDiffMarkers_combined[[8]],listDiffMarkers_combined[[9]])


## Save results ##
saveRDS(listDiffMarkers_combined,file=paste0(sampleFolder,"Robjects/markersDiffSamples_combined_9pops.rds"))
saveRDS(listDiffMarkers_all,file=paste0(sampleFolder,"Robjects/markersDiffSamples_all_in_one_9pops.rds"))

### Write to Excel
library('openxlsx')
write.xlsx(listDiffMarkers_combined, file = paste0(sampleFolder,"results/summaryDiffMarkers_combined_9pops.xlsx"))
write.xlsx(listDiffMarkers_all, file = paste0(sampleFolder,"results/summaryDiffMarkers_all_in_one_9pops.xlsx"))

## Read results ##
listDiffMarkers_combined<-readRDS(file=paste0(sampleFolder,"Robjects/markersDiffSamples_combined_9pops.rds"))
listDiffMarkers_all<-readRDS(file=paste0(sampleFolder,"Robjects/markersDiffSamples_all_in_one_9pops.rds"))

# Remove DE genes from Fras1+ LV from figure! Only 19 cells LpsNeg and 11 cells LpsPos! Not trustworthy!
listDiffMarkers_combined$Xist_CPE<-listDiffMarkers_combined$Xist_CPE[-which(listDiffMarkers_combined$Xist_CPE$Location == "LV"),]

#Reorder lists: 
listDiffMarkers_combined<-listDiffMarkers_combined[c(1,2,4,7,5,6,3,9,8)]

Titles<-as.factor(c('Epithelial Cells', "Fras1+ Epithelial Cells",
                    "Endothelial Cells", "Vascular Assocd. Cells", "Fibroblasts Type 1", "Fibroblasts Type 2", 
                    "Macrophages", "Dendritic Cells", "NK Cells"))

Cell_pops_plot<-c("Epithelial Cells","Xist+ Epithelial Cells","Endothelial Cells","Vascular Associated Cells",
                  "Fibroblasts Type 1", "Fibroblasts Type 2","Macrophages","Dendritic Cells","NK Cells") #Order Roos

#New plot Daan
Colorset_blind<-colorblind_pal()(8)
Colorset_blind[8]<-"#ff08e8"

## Combined expression and gene class ##
myList=list()
myList2<-list()

for (i in 1:length(listDiffMarkers_combined)) {
  testData<-listDiffMarkers_combined[[i]]
  testData<-testData[,c(2,6,7,8,14)]
  testData[,"Jitter"]<-runif(length(testData$geneSymbol)) #Add random numbers for x axis uniform distribution (Jitter to see dots)
  testData[which(testData$Location=="4V"),"Jitter"]<-testData[which(testData$Location=="4V"),"Jitter"]+1.5
  testData1<-testData[which(testData$Location=="LV"),]
  testData2<-testData[which(testData$Location=="4V"),]
  #Expression data table
  testData1[,"Expression"]<- apply(seuratObjNew2@assays$RNA@data[testData1[,2],rownames(seuratObjNew2@meta.data[which(Idents(seuratObjNew2)==paste0(Cell_pops_plot[i],"_LpsPosLat")),])],1,function(x){mean(x)}) #Apply mean over the rows (for each gene over all cells for that cell population)
  testData2[,"Expression"]<- apply(seuratObjNew2@assays$RNA@data[testData2[,2],rownames(seuratObjNew2@meta.data[which(Idents(seuratObjNew2)==paste0(Cell_pops_plot[i],"_LpsPosFour")),])],1,function(x){mean(x)}) #Apply mean over the rows (for each gene over all cells for that cell population)
  testData<-rbind(testData1,testData2)
  testData<-testData[order(testData$Expression),] #Order on expression
  for (k in 1:nrow(testData)){
    if(testData[k,"avg_logFC"]<0){
      testData[k,"Expression"]<-testData[k,"Expression"]*-1
    }
  }
  #Gene class data table
  testDataClass<-testData
  testDataClass<-testDataClass[order(testDataClass$Gene_class),] #Order on gene class

  if (i==1){
    # Expression plot
    p <- ggplot(testData, aes(x=Jitter, y=avg_logFC, label = geneSymbol, color=Expression)) + 
      geom_point(shape=16) + 
      scale_color_gradient2(midpoint=0, low="blue", mid="white",
                            high="red", space ="Lab" ) +
      coord_cartesian(ylim = c(-3, 6)) +
      geom_text(aes(x=2.2, label="4V", y=6.3), color="lightblue", hjust=0) +
      geom_text(aes(x=-0.2, label="LV", y=6.3), color="lightblue", hjust=0) +
      geom_vline(xintercept = 1.25, color = "lightblue", linetype = "longdash", size = 1) +
      geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-0.24,ymax=0.24),alpha=0.01,color = NA)+
      labs(title=as.character(Titles[i])) +
      
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size=14)) #,legend.position = "none"
    
    # Gene class plot
    testDataClass$Gene_class<-as.factor(testDataClass$Gene_class)
    levels(testDataClass$Gene_class)<-Colorset_blind[c(1,2,3,8)] #Choose 4 levels!
    
    p2 <- ggplot(testDataClass, aes(x=Jitter, y=avg_logFC, label = geneSymbol)) + 
      geom_point(shape=21, fill=testDataClass$Gene_class) + 
      coord_cartesian(ylim = c(-3, 6)) +
      geom_text(aes(x=2.2, label="4V", y=6.3), color="lightblue", hjust=0) +
      geom_text(aes(x=-0.2, label="LV", y=6.3), color="lightblue", hjust=0) +
      geom_vline(xintercept = 1.25, color = "lightblue", linetype = "longdash", size = 1) +
      geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-0.24,ymax=0.24),alpha=0.01,color = NA)+
      labs(title=as.character(Titles[i])) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size=14))
    
    
  } else if (i != 2 && i !=6) { #Have all 4 levels ##(i != 5 && i !=7) #Order Daan (i != 2 && i !=6) #Order Roos
    # Expression plot
    p <- ggplot(testData, aes(x=Jitter, y=avg_logFC, label = geneSymbol, color=Expression)) + 
      geom_point(shape=16) +
      scale_color_gradient2(midpoint=0, low="blue", mid="white",
                            high="red", space ="Lab" ) +
      coord_cartesian(ylim = c(-3, 6)) +
      geom_text(aes(x=2.2, label="4V", y=6.3), color="lightblue", hjust=0) +
      geom_text(aes(x=-0.2, label="LV", y=6.3), color="lightblue", hjust=0) +
      geom_vline(xintercept = 1.25, color = "lightblue", linetype = "longdash", size = 1) +
      geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-0.24,ymax=0.24),alpha=0.01,color = NA)+
      labs(title=as.character(Titles[i])) +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size=14))
    
    # Gene class plot
    testDataClass$Gene_class<-as.factor(testDataClass$Gene_class)
    levels(testDataClass$Gene_class)<-Colorset_blind[c(1,2,3,8)] #Choose 4 levels!
    
    p2 <- ggplot(testDataClass, aes(x=Jitter, y=avg_logFC, label = geneSymbol)) + 
      geom_point(shape=21, fill=testDataClass$Gene_class) + 
      coord_cartesian(ylim = c(-3, 6)) +
      geom_text(aes(x=2.2, label="4V", y=6.3), color="lightblue", hjust=0) +
      geom_text(aes(x=-0.2, label="LV", y=6.3), color="lightblue", hjust=0) +
      geom_vline(xintercept = 1.25, color = "lightblue", linetype = "longdash", size = 1) +
      geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-0.24,ymax=0.24),alpha=0.01,color = NA)+
      labs(title=as.character(Titles[i])) +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size=14))
    
  } else { #Only have 3 levels
    # Expression plot
    p <- ggplot(testData, aes(x=Jitter, y=avg_logFC, label = geneSymbol, color=Expression)) + 
      geom_point(shape=16) + 
      scale_color_gradient2(midpoint=0, low="blue", mid="white",
                            high="red", space ="Lab" ) +
      coord_cartesian(ylim = c(-3, 6), xlim = c(-0.2,2.5)) +
      geom_text(aes(x=2.2, label="4V", y=6.3), color="lightblue", hjust=0) +
      geom_text(aes(x=-0.2, label="LV", y=6.3), color="lightblue", hjust=0) +
      geom_vline(xintercept = 1.25, color = "lightblue", linetype = "longdash", size = 1) +
      geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-0.24,ymax=0.24),alpha=0.01,color = NA)+
      labs(title=as.character(Titles[i])) +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size=14)) 
    
    # Gene class plot
    testDataClass$Gene_class<-as.factor(testDataClass$Gene_class)
    levels(testDataClass$Gene_class)<-Colorset_blind[c(1,2,8)] #Choose 3 levels
    
    p2 <- ggplot(testDataClass, aes(x=Jitter, y=avg_logFC, label = geneSymbol)) + 
      geom_point(shape=21, fill=testDataClass$Gene_class) + 
      coord_cartesian(ylim = c(-3, 6), xlim = c(-0.2,2.5)) +
      geom_text(aes(x=2.2, label="4V", y=6.3), color="lightblue", hjust=0) +
      geom_text(aes(x=-0.2, label="LV", y=6.3), color="lightblue", hjust=0) +
      geom_vline(xintercept = 1.25, color = "lightblue", linetype = "longdash", size = 1) +
      geom_rect(aes(xmin=-Inf,xmax=Inf,ymin=-0.24,ymax=0.24),alpha=0.01,color = NA)+
      labs(title=as.character(Titles[i])) +
      theme(axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5,size=14)) 
  }
  myList[[i]]<-p
  myList2[[i]]<-p2
}


ggsave(grid.arrange(grobs = myList, ncol = length(Titles), nrow = 1), file=paste0(sampleFolder,"results/QC/17_DE_gene_plot_final_Expression.png"), height= 8, width = 30, dpi = 300)
ggsave(grid.arrange(grobs = myList2, ncol = length(Titles), nrow = 1), file=paste0(sampleFolder,"results/QC/17_DE_gene_plot_final_Gene_Class.png"), height= 8, width = 20, dpi = 300)


#####################################################################################################################

## GO enrichment analysis with sc background 
Background_sc<-rownames(seuratObjNew2) 
universe<-Background_sc

## Check GO enrichment within BBCs (LpsPos vs LpsNeg) 
## Look at diff markers 4V and LV combined 
listDiffMarkers<-readRDS(file=paste0(sampleFolder,"Robjects/markersDiffSamples_Full.rds"))

## BBCs
All_DE_in_Stalk<-enrichGO(
  as.character(listDiffMarkers$Fibroblasts_Type_2_LpsPosvsLpsNeg$geneSymbol),
  'org.Mm.eg.db',
  keyType = "SYMBOL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  universe,
  qvalueCutoff = 0.2,
  minGSSize = 10,
  maxGSSize = 500,
  readable = FALSE,
  pool = T
)

## Create dotplots
D3<-dotplot(All_DE_in_Stalk, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")

## Save dotplots
pdf(file=paste0(sampleFolder,"results/QC/Stalk_FB_analysis/EnrichGO_Dotplot_ALL_SCbackground_",sampleName,".pdf"), width = 10, height = 10)
D3
dev.off()

## Save results
saveRDS(All_DE_in_Stalk, paste0(sampleFolder,"results/QC/Stalk_FB_analysis/EnrichGO_Dotplot_ALL_SCbackground_",sampleName,".rds"))
write.xlsx(All_DE_in_Stalk@result,file=paste0(sampleFolder,"results/QC/Stalk_FB_analysis/EnrichGO_results_ALL_SCbackground_",sampleName,".xlsx"))

## Filter to only BP
All_DE_in_Stalk_filtered<-All_DE_in_Stalk
All_DE_in_Stalk_filtered@result<-All_DE_in_Stalk_filtered@result[-which(All_DE_in_Stalk_filtered@result$ONTOLOGY == "MF" | All_DE_in_Stalk_filtered@result$ONTOLOGY == "CC"),]
D3_filtered<-dotplot(All_DE_in_Stalk_filtered, split="ONTOLOGY", showCategory = 30 ) + facet_grid(ONTOLOGY~., scale="free")
pdf(file=paste0(sampleFolder,"results/QC/Stalk_FB_analysis/EnrichGO_Dotplot_ALL_SCbackground_only_BP_",sampleName,".pdf"), width = 15, height = 10)
D3_filtered
dev.off()
write.xlsx(All_DE_in_Stalk_filtered@result,file=paste0(sampleFolder,"results/QC/Stalk_FB_analysis/EnrichGO_results_ALL_SCbackground_only_BP_",sampleName,".xlsx"))

#########################################################################################################

## Check gene list chemokines 
Daan_list<-read.xlsx("1.Documentation/Chemokine list for ChP BBCs.xlsx", sheet = 2)

###First letter upper case
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

Figure_genes<-firstup(tolower(Daan_list$Chemokines))
Figure_genes<-intersect(Figure_genes,rownames(seuratObjNew2))

# Create dotplot, featureplots and violoin plots
Colors_dotplot<-c("#071AE5","#F50635") #030720
D1<-DotPlot(seuratObjNew, features = rev(Figure_genes), cols = Colors_dotplot) + RotatedAxis()

pdf(file=paste0(sampleFolder,"results/QC/Feature_plots_and_dotplots_chemokines_2023_",sampleName,".pdf"), height = 10, width = 15)
D1

for (feature in Figure_genes) {
  F1<-FeaturePlot(object = seuratObjNew, features =feature, cols = c("yellow", "red"), 
                  reduction = "umap", min.cutoff = 'q2', max.cutoff = 'q98', label = T, repel = T, pt.size = 1.5, order=T)
  print(F1)
}

dev.off()

pdf(file=paste0(sampleFolder,"results/QC/Violinplots_chemokines_2023_",sampleName,".pdf"), height = 10, width = 20)

for (feature in c(Figure_genes,"Dpp4")) {
  F1<-  VlnPlot(seuratObjNew2, features = feature, col=Colorset)  + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(F1)
  F2<-  VlnPlot(seuratObjNew2, features = feature, col=Colorset[c(11,12)], idents = levels(Idents(seuratObjNew2))[c(11,12)])  + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(F2)
}

dev.off()
