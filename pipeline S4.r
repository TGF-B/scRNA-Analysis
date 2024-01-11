
#now you have the barcodes.tsv, features.tsv, matrix.mtx after running the cellranger

####################################build/activate the virtual environment###################

#conda create -n myenv r-base
conda activate myenv
#enter R
R
#virtual environment keep failing, so excuted locally
###########################################load the data##########################################
#import the the extracted datas and merge
library(Seurat) #version 4.4.0
library(dplyr) #helps to rename the column
library(Matrix) #Matrix 1.6.1 or there will be an as_cholmod_sparse error
library(patchwork) #helps to plot multiple plots in one figure
#change the loading command as you want
#  ./xxx  consits of three files: barcodes.tsv, features.tsv, matrix.mtx
CP1.data<-Read10X(data.dir="./CP1")
CP2.data<-Read10X(data.dir="./CP2")
CP3.data<-Read10X(data.dir="./CP3")
Ctrl1.data<-Read10X(data.dir="./Ctrl1")
Ctrl2.data<-Read10X(data.dir="./Ctrl2")
Ctrl3.data<-Read10X(data.dir="./Ctrl3")
group<-list(CP1=CP1.data,CP2=CP2.data,CP3=CP3.data,Ctrl1=Ctrl1.data,Ctrl2=Ctrl2.data,Ctrl3=Ctrl3.data)
names(group)<-c("CP1","CP2","CP3","Ctrl1","Ctrl2","Ctrl3")
save(group,file="./RData/group.RData")

Filter<-function(CP1){
CP1<-CreateSeuratObject(counts=CP1.data,project='MouseCP1',min.cells=3,min.features=200)
CP1[["percent_mt"]]<-PercentageFeatureSet(CP1,pattern="^mt-")
#refine the data
CP1<-subset(CP1,subset=nFeature_RNA>200&nFeature_RNA<2500&percent_mt<5)
VlnPlot(CP1,features=c("nFeature_RNA","nCount_RNA","percent_mt"),ncol=3)
#pca
CP1<-NormalizeData(CP1)
CP1<-FindVariableFeatures(CP1,selection.method="vst",nfeatures=2000)
all.genes<-rownames(CP1)
CP1<-ScaleData(CP1,features=all.genes)
plot1<-RunPCA(CP1,features=VariableFeatures(object=CP1))
plot2<-FeatureScatter(CP1,feature1="nCount_RNA",feature2="nFeature_RNA")
plot3<-FeatureScatter(CP1,feature1="nCount_RNA",feature2="percent_mt")
#show the 3 plots
plot1+plot2+plot3
return(CP1)
}

#laply filter function to the whole group list
loaded_group<-lapply(group,FUN=Filter)
names(loaded_group)<-c("CP1","CP2","CP3","Ctrl1","Ctrl2","Ctrl3")
#check the raw data quality

'
Rename the elements in Orig.ident column as CPX
CP1@meta.data <- CP1@meta.data %>%
 dplyr::mutate(orig.ident = ifelse(orig.ident == "SeuratProject", "CP1", orig.ident))
CP2@meta.data <- CP2@meta.data %>%
 dplyr::mutate(orig.ident = ifelse(orig.ident == "SeuratProject", "CP2", orig.ident))
CP3@meta.data <- CP3@meta.data %>%
 dplyr::mutate(orig.ident = ifelse(orig.ident == "SeuratProject", "CP3", orig.ident))
Ctrl1@meta.data <- Ctrl1@meta.data %>%
 dplyr::mutate(orig.ident = ifelse(orig.ident == "SeuratProject", "Ctrl1", orig.ident))
Ctrl2@meta.data <- Ctrl2@meta.data %>%
 dplyr::mutate(orig.ident = ifelse(orig.ident == "SeuratProject", "Ctrl2", orig.ident))
Ctrl3@meta.data <- Ctrl3@meta.data %>%
 dplyr::mutate(orig.ident = ifelse(orig.ident == "SeuratProject", "Ctrl3", orig.ident))
'
save(loaded_group,file="./RData/loaded_group.RData")

#############pre process#############################
dims<-15 #decided from ElbowPlot(omitted)

library(DoubletFinder)#version 2.0.3 
#paramSweep_v3.R and DoubletFinder_v3.R have been modified ,where meta@data$RNA@counts is changed to GetAssayData(meta, "counts")
library(tidyverse)
library(patchwork)

Pre_process<-function(data){
data<-NormalizeData(object=data)
data<-FindVariableFeatures(data,selection.method="vst",nfeatures=2000)
all.genes<-rownames(data)
data<-ScaleData(data,features=all.genes)
data<-RunPCA(data,features=VariableFeatures(object=data))
VizDimLoadings(data,reduction="pca",dims=1:dims)   
DimHeatmap(data,dims=1:dims,cells=500,balanced=TRUE)
data<-FindNeighbors(data,dims=1:dims)
data<-FindClusters(data,resolution=0.6)
data<-RunTSNE(data,dims=1:dims)
#data<-RunUMAP(data,dims=1:15) #cannont run UMAP because as_cholmod_sparse not provided by package "Matrix"
	return(data)
}
processed_group<-lapply(loaded_group,FUN=Pre_process) # processe the  data and save them in a list
names(processed_group) <- paste0("processed_",names(processed_group)) #add the prefix "processed_" to the names of the list
save(processed_group,file="RData/processed_group.RData")

##########################de doublet cells###############################

Find_doublet <- function(processed_xxx) {
sweep.res.list <- paramSweep_v3(processed_xxx, PCs=1:dims, sct = FALSE, num.cores=1)#seu@assays$RNA@counts to Seurat::GetAssayData(seu, "counts")
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #choose the best pK
DoubletRate = ncol(processed_xxx) * 8 * 1e-6 #0.027296
nExp_poi <- round(DoubletRate * ncol(processed_xxx)) #93
#homotypic.prop <- modelHomotypic(data, PCs = 1:dims, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
processed_xxx <- doubletFinder_v3(processed_xxx, PCs = 1:dims, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#colnames(processed_xxx@meta.data)[colnames(processed_xxx@meta.data) == "pANN_0.25_0.05_93"] <- "Double_score"
#colnames(processed_xxx@meta.data)[colnames(processed_xxx@meta.data) == "DF.classifications_0.25_0.05_93"] <- "doublet_info" 
#pANN_0.25_0.05_93(doublet score) and DF.classifications_0.25_0.05_93(doublet or singlet) are column names in data matrix, 
#stands for the parameter pN=0.25,pK=0.05,nExp=93, the name depends on the parameter of each run
######## So it would be better to choose columne by its column number instead of its name #######
colnames(processed_xxx@meta.data)[ncol(processed_xxx@meta.data)] = "doublet_info" #ncol returns column number, the DF.xxx is the last column
DimPlot(processed_xxx,reduction="tsne",label=T,group.by="doublet_info")
processed_xxx<- subset(processed_xxx, subset = doublet_info == "Singlet")
return(processed_xxx)
}

de2_group<-lapply(processed_group,FUN=Find_doublet)
names(de2_group) <- c("de2_CP1","de2_CP2","de2_CP3","de2_Ctrl1","de2_Ctrl2","de2_Ctrl3")

save(de2_group,file="./RData/de2_group.RData")

#################################################Merge###########################
#merge the three de2_group
de2_group$de2_CP1@meta.data$orig.ident <- "CP1"
de2_group$de2_CP2@meta.data$orig.ident <- "CP2"
de2_group$de2_CP3@meta.data$orig.ident <- "CP3"
de2_group$de2_Ctrl1@meta.data$orig.ident <- "Ctrl1"
de2_group$de2_Ctrl2@meta.data$orig.ident <- "Ctrl2"
de2_group$de2_Ctrl3@meta.data$orig.ident <- "Ctrl3"

Merged_CP<- merge(de2_group$de2_CP1, y = list(de2_group$de2_CP2, de2_group$de2_CP3))
Merged_Ctrl<- merge(de2_group$de2_Ctrl1, y = list(de2_group$de2_Ctrl2, de2_group$de2_Ctrl3))

save(Merged_CP,file="RData/Merged_CP.RData")
save(Merged_Ctrl,file="RData/Merged_Ctrl.RData")

##############################de batch effect####################################
library(harmony) #version 1.2.0

debatch<-function(Merged_xxx){
Merged_xxx<-NormalizeData(Merged_xxx)
Merged_xxx<-FindVariableFeatures(Merged_xxx)
all.genes_xxx<-rownames(Merged_xxx)
Merged_xxx<-ScaleData(Merged_xxx,features=all.genes_xxx)
Merged_xxx<-RunPCA(Merged_xxx,npcs=50,verbose=FALSE)
Merged_xxx_a<-RunHarmony(Merged_xxx,"orig.ident")
Merged_xxx<-RunUMAP(Merged_xxx_a,dims=1:15,reduction="harmony")
#DimPlot(Merged_xxx,reduction="umap",label=T,group.by="orig.ident")
Merged_xxx<-FindNeighbors(Merged_xxx,dims=1:15,reduction="harmony")
Merged_xxx<-FindClusters(Merged_xxx,resolution=0.4,reduction="harmony")
DimPlot(Merged_xxx,reduction="umap",label=T,group.by="orig.ident")
return(Merged_xxx)
}
Merged_CP_debatched<-debatch(Merged_CP)
Merged_Ctrl_debatched<-debatch(Merged_Ctrl)

save(Merged_CP_debatched,file="RData/Merged_CP_bebatched.RData")
save(Merged_Ctrl_debatched,file="RData/Merged_Ctrl_debatched.RData")

###########################Quality Control#####################################

QC<-function(Merged_xxx){
#mt RNA percentage
plot1 <- FeatureScatter(Merged_xxx, feature1 = "nCount_RNA", feature2 = "percent_mt")
plot2 <- FeatureScatter(Merged_xxx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#High variable features
Merged_xxx <- FindVariableFeatures(Merged_xxx,selection.method = "vst", nfeatures = 2000, verbose = FALSE)
plot3 <- VariableFeaturePlot(Merged_xxx)
top30 <-head(VariableFeatures(Merged_xxx),30)
plot3 <- LabelPoints(plot = plot3, points = top30, repel = TRUE, xnudge = 0, ynudge = 0)
#scale the data(ignore the overlap error when occurs)
all.genes <- rownames(Merged_xxx)
Merged_xxx <- ScaleData(Merged_xxx, features = all.genes)
Merged_xxx<-RunPCA(Merged_xxx,features=VariableFeatures(object = Merged_xxx))
plot4<-DimPlot(Merged_xxx, reduction = "pca", label = TRUE, pt.size = 0.5) + NoLegend()
# Examine and visualize PCA results in a few different ways
#print(Merged_xxx[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Merged_xxx, dims = 1:15, reduction = "pca")
plot5<-DimPlot(Merged_xxx, reduction = "pca") + NoLegend()
#Cluster the cells
Clustering <- FindNeighbors(Merged_xxx, dims = 1:15)
Clustered <- FindClusters(Clustering, resolution = 0.5)
Merged_xxx.markers <- FindAllMarkers(Merged_xxx, only.pos = TRUE)
Merged_xxx.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1)
#DimHeatPlot
top30 <- head(VariableFeatures(Merged_xxx), 10)
Merged_xxx.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 5) %>% #5 percent of all marker genes
    ungroup() %>%
    dplyr::filter(gene %in% top30)
DoHeatmap(Clustered, features = top30) + NoLegend()
return(Merged_xxx)
}

Merged_CP_QC<-QC(Merged_CP_debatched)
Merged_Ctrl_QC<-QC(Merged_Ctrl_debatched)
save(QC_Merged_CP,file="RData/Merged_CP_QC.RData")
save(QC_Merged_Ctrl,file="RData/Merged_Ctrl_QC.RData")

###############################automated Annotation with SingleR#################################

library(SingleR) #version 2.2.0
library(celldex) #version 1.10.1
library(pheatmap)
library(ggplot2)

#turning Seurat matrix into SingleR matrix
SingleR_xxx<-GetAssayData(Merged_xxx,slot="data")
#import mouse gene marker library manually(cuz automated failed)
MouseRNAseq<-readRDS("MouseRNAseqData.rds") 
#general match to one dataset
Merged_xxx.match<-SingleR(test=SingleR_xxx,ref=MouseRNAseq,labels=MouseRNAseq$label.main)
#SingleR table PCA and define the labels
table(Merged_xxx.match$labels,Merged_xxx$seurat_clusters)
Merged_xxx@meta.data$labels <-Merged_xxx.match$labels
#PCA
Merged_PCA<-DimPlot(Merged_xxx,group.by=c("labels"),reduction="pca")
#dims=c(2,1),rotate the plot because unexplainale data appearance
#remove the duplicates
Merged_xxx<-subset(Merged_xxx,subset=labels!="NA")
#TSNE
Merged_TSNE<- RunTSNE(Merged_xxx, dims = 1:15)

plot_CP<-DimPlot(Merged_TSNE, group.by=c("labels"),reduction = "tsne",dims=c(2,1))

plot_Ctrl<-DimPlot(Merged_TSNE, group.by=c("labels"),reduction = "tsne",dims=c(2,1))
#show the DimPlot

#manually Annatation

#Multiple Violin pay attention to the dims match for new.cluster and Merged_Ctrl
#new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
#names(new.cluster.ids) <- levels(Merged_Ctrl)
#Merged_Ctrl <- RenameIdents(Merged_Ctrl, new.cluster.ids)
#DimPlot(Merged_Ctrl, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#library(ggplot2)
#plot <- DimPlot(Merged_Ctrl_Clustered, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
#    theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
#ggsave(filename = "./Ctrl_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
#saveRDS(Merged_Ctrl, file = "./")

###########################################plot the de-dims############################

#UMAP unbale to comply because as_cholmod_sparse not provided by package "Matrix")
#Merged_CP_Clustered_UMAP <- RunUMAP(Merged_CP_Clustered, dims = 1:13)
#plot_CP_UMAP<-DimPlot(Merged_CP_Clustered_UMAP, group.by=c("labels"),reduction = "umap")

#Merged_Ctrl_Clustered_UMAP <- RunUMAP(Merged_Ctrl_Clustered, dims = 1:13)
#plot_Ctrl_UMAP<-DimPlot(Merged_Ctrl_Clustered_UMAP, group.by=c("Ctrl"),reduction = "tsne",dims=c(2,1))
#plot_CP_UMAP + plot_Ctrl_UMAP
