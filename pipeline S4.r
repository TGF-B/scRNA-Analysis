
#now you have the barcodes.tsv, features.tsv, matrix.mtx after running the cellranger

######remember using "str" often to check the structure whenever come across a new data #######

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
CP1<-Read10X(data.dir="./CP1")
CP2<-Read10X(data.dir="./CP2")
CP3<-Read10X(data.dir="./CP3")
Ctrl1<-Read10X(data.dir="./Ctrl1")
Ctrl2<-Read10X(data.dir="./Ctrl2")
Ctrl3<-Read10X(data.dir="./Ctrl3")
group<-list(CP1,CP2,CP3,Ctrl1,Ctrl2,Ctrl3)
names(group)<-c("CP1","CP2","CP3","Ctrl1","Ctrl2","Ctrl3")
save(group,file="./RData/group.RData")
#load("./RData/group.RData")
Filter<-function(xxx){
xxx<-CreateSeuratObject(counts=xxx,project='MouseCP1',min.cells=3,min.features=200)
xxx[["percent_mt"]]<-PercentageFeatureSet(xxx,pattern="^mt-")
#refine the data
xxx<-subset(xxx,subset=nFeature_RNA>200&nFeature_RNA<2500&percent_mt<5)
VlnPlot(xxx,features=c("nFeature_RNA","nCount_RNA","percent_mt"),ncol=3)
#pca
xxx<-NormalizeData(xxx)
xxx<-FindVariableFeatures(xxx,selection.method="vst",nfeatures=2000)
all.genes<-rownames(xxx)
xxx<-ScaleData(xxx,features=all.genes)
RunPCA(xxx,features=VariableFeatures(object=xxx))
FeatureScatter(xxx,feature1="nCount_RNA",feature2="nFeature_RNA")
FeatureScatter(xxx,feature1="nCount_RNA",feature2="percent_mt")
return(xxx)
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
load("RData/loaded_group.RData")
dims<-8 #decided from ElbowPlot(omitted)

library(DoubletFinder)#version 2.0.3 
#paramSweep_v3.R and DoubletFinder_v3.R have been modified ,where meta@data$RNA@counts is changed to GetAssayData(meta, "counts")
#try devtools::install("./DoubletFinder-master") locally
#library(tidyverse)
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
data<-RunUMAP(data,dims=1:dims)
#plot the UMAP
DimPlot(data,reduction="umap"),label=T,group.by="orig.ident" )
#data<-RunUMAP(data,dims=1:15) #sometimes cannont run UMAP because as_cholmod_sparse not provided by package "Matrix"
	return(data)
}
processed_group<-lapply(loaded_group,FUN=Pre_process) # processe the  data and save them in a list
#names(processed_group) <- paste0("processed_",names(processed_group)) #add the prefix "processed_" to the names of the list
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
save(de2_group,file="./RData/de2_group.RData")
names(de2_group) <- c("de2_CP1","de2_CP2","de2_CP3","de2_Ctrl1","de2_Ctrl2","de2_Ctrl3")


#################################################Merge###########################
#merge the three de2_group
load("RData/de2_group.RData")
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
load("RData/Merged_CP.RData")
load("RData/Merged_Ctrl.RData")
#Merged_CP=subset(Merged_CP,subset=orig.ident!="SeuratProject")
#Merged_Ctrl=subset(Merged_Ctrl,subset=orig.ident!="SeuratProject")
library(harmony) #version 1.2.0
#devtools::install_github("immunogenomics/harmony", ref = "develop")

debatch<-function(Merged_xxx){
Merged_xxx<-NormalizeData(Merged_xxx)
Merged_xxx<-FindVariableFeatures(Merged_xxx)
all.genes_xxx<-rownames(Merged_xxx)
Merged_xxx<-ScaleData(Merged_xxx,features=all.genes_xxx)
Merged_xxx<-RunPCA(Merged_xxx,npcs=50,verbose=FALSE)
Merged_xxx<-RunHarmony(Merged_xxx,"orig.ident")
Merged_xxx<-RunUMAP(Merged_xxx,dims=1:dims,reduction="harmony")#error because of Matrix package(as_cholmod_sparse)
#Merged_xxx<-RunTSNE(Merged_xxx,dims=1:15,reduction="harmony")
#DimPlot(Merged_xxx,reduction="umap",label=T,group.by="orig.ident")
Merged_xxx<-FindNeighbors(Merged_xxx,dims=1:dims,reduction="harmony")
Merged_xxx<-FindClusters(Merged_xxx,resolution=0.4,reduction="harmony")
DimPlot(Merged_xxx,reduction="umap",label=T,group.by="orig.ident")
return(Merged_xxx)
}

Merged_Ctrl_debatched<-debatch(Merged_Ctrl)
plot7<-DimPlot(Merged_Ctrl_debatched,reduction='umap',label=T,group.by='orig.ident')
#plot7<-DimPlot(Merged_Ctrl_debatched,reduction='tsne',label=T,group.by='orig.ident')
Merged_CP_debatched<-debatch(Merged_CP)
plot8<-DimPlot(Merged_CP_debatched,reduction='umap',label=T,group.by='orig.ident')
#plot8<-DimPlot(Merged_CP_debatched,reduction='tsne',label=T,group.by='orig.ident')
plot7+plot8

save(Merged_CP_debatched,file="RData/Merged_CP_debatched.RData")
save(Merged_Ctrl_debatched,file="RData/Merged_Ctrl_debatched.RData")

###########################Quality Control(Optional)#####################################
load("RData/Merged_Ctrl_debatched.RData")
load("RData/Merged_CP_debatched.RData")
#High variable features
Merged_xxx<-Merged_Ctrl_debatched
#scatter the variable features
FeatureScatter(Merged_xxx,feature1="nCount_RNA",feature2="nFeature_RNA")
#FeatureScatter(Merged_xxx,feature1="nCount_RNA",feature2="percent.mt")
Merged_xxx <- FindVariableFeatures(Merged_xxx,selection.method = "vst", nfeatures = 2000, verbose = FALSE)
plot9 <- VariableFeaturePlot(Merged_xxx)
top30 <-head(VariableFeatures(Merged_xxx),30)
plot9 <- LabelPoints(plot9, points = top30, repel = TRUE, xnudge = 0, ynudge = 0)

Merged_xxx2<-Merged_CP_debatched
plot<-FeatureScatter(Merged_xxx,feature1="nCount_RNA",feature2="nFeature_RNA")
#plot<-FeatureScatter(Merged_xxx,feature1="nCount_RNA",feature2="percent.mt")
Merged_xxx2<- FindVariableFeatures(Merged_xxx2,selection.method = "vst", nfeatures = 2000, verbose = FALSE)
plot10 <- VariableFeaturePlot(Merged_xxx2)
top30_2 <-head(VariableFeatures(Merged_xxx2),30)
plot10 <- LabelPoints(plot10, points = top30_2, repel = TRUE, xnudge = 0, ynudge = 0)

plot9+plot10 #the same

#scale the data(ignore the overlap error when occurs)
top30 <- head(VariableFeatures(Merged_xxx), 10)
Merged_xxx.markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 5) %>% #5 percent of all marker genes
    ungroup() %>%
    dplyr::filter(gene %in% top30)
DoHeatmap(Clustered, features = top30) + NoLegend()
save(CP_QC,file="RData/Merged_CP_QC.RData")
save(Ctrl_QC,file="RData/Merged_Ctrl_QC.RData")

#cluster by markers manually.........
Merged_debatched<-merge(Merged_Ctrl_debatched,Merged_CP_debatched)
Pre_process(Merged_debatched)
Merged_debatched<-debatch(Merged_debatched)
'
features_Monocytes<-c("Apoe","Cd14","Cd68","C5ar1","Csf1r","Lyz2","Ly6c2","Ccr2","Mgst1")#cell IV
RidgePlot(Merged_debatched,features=features_Monocytes,ncol=3)
features_NK_cells<-c("Ccl5","Fos","Gnly","Nkg7","Ccl4","Il7r")#cell tracking
NK_detect<-RidgePlot(Merged_debatched,features=features_NK_cells,ncol=3)
ggsave("NK_detect.pdf",plot=NK_detect,width=10,height=10)
features_DC_cells<-c("Xcr1","clec9a","Ccl22","Ccr7","Basp1","Bst2","Cd11c","Cd16","Ccl5","ly6d","S100a4","Tcf7")#cell tracking
DC_detect<-RidgePlot(Merged_debatched,features=features_DC_cells,ncol=3)
ggsave("DC_detect.pdf",plot=DC_detect,width=10,height=10)
'

features_new<-c("Ncr1","Dnajb1","Crem","Bcl2l11","Nkg7","Prf1","Icos","Zbtb16","Cntn1","Cd3g","Cd68","Chil3","Ccr2","Itgax","Batf3","Bst2","Cd19","Stat1","Cd31","Prox1","Cd25","Ccl19")#from lymph node(2 articles with Cellmarker)
cluster_marker<-DotPlot(Merged_debatched,features=features_new)+RotatedAxis()
ggsave("cluster_marker.pdf",plot=cluster_marker,width=10,height=10)
DimPlot(Merged_debatched,reduction='umap',label=T)

Merged<-Merged_debatched
#define the cluster with celltype manually
new.cluster.id<-c(
  "NK cells","Dendritic cells","NK cells","NK cells","Dendritic cells","T cells","NK cells","NK cells","Mono/Mac",
  "Dendritic cells","B cells","Endothelial cells","T cells","Fibroblasts"
)
names(new.cluster.id)<-levels(Merged)
Merged_new<-RenameIdents(Merged_debatched,new.cluster.id)
Merged_debatched$celltype<-Idents(Merged_new)

cluster_color_manually<-c("NK cells"="#836621","T cells"="#c9c941","Mono/Mac"="#ae4ef3","Dendritic cells"="#5ac8c9",
"B cells"="#7a73f4","Endothelial cells"="#5cc86e","Fibroblasts"="#e83123")
Merged_debatched_data<-Merged_debatched@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=Merged_debatched$celltype)
axis<-Merged_debatched_data %>%
  group_by(celltype) %>%
   dplyr::summarise(
      UMAP_1=median(UMAP_1),
      UMAP_2=median(UMAP_2))

Manually_UMAP<-DimPlot(Merged_debatched, reduction = "umap",group.by="celltype",pt.size=0.3,label=F)+scale_color_manual(values=cluster_color_manually)+
theme_dr(xlength=0.3,ylength=0.3,arrow=grid::arrow(length=unit(0.2,"inches"),type="closed"))+theme(panel.grid=element_blank())+geom_label_repel(data=axis,aes(x=UMAP_1,y=UMAP_2,label=celltype),fontface="bold",point.padding=unit(0.3,"lines"))
ggsave("Manually_UMAP.pdf",plot=Manually_UMAP,width=10,height=10)

cluster_marker<-DotPlot(Merged_debatched,features=features_new)+RotatedAxis()
ggsave("cluster_marker.pdf",plot=cluster_marker,width=10,height=10)
DimPlot(Merged_debatched,reduction='umap',label=T)
###############################automated Annotation with SingleR#################################
#load RData/Merged_CP_QC.RData

#####debatched part wrong#####
#####debatched part wrong#####
#####debatched part wrong#####

library(SingleR) #version 2.4.1 depended on Matrix 1.6-1
#install.packages("SingleR_1.4.1.tar.gz",repos = NULL,type="source")
library(celldex) #version 1.10.1
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

library(ggsci)

#import mouse gene marker library manually(cuz automated always failed)
MouseRNAseq<-readRDS("MouseRNAseqData.rds") 
Ctrl<-Merged_Ctrl_debatched
CP<-Merged_CP_debatched

######plot before annotation
#plot11<-DimPlot(Ctrl, reduction = "tsne",label=T)
#plot12<-DimPlot(CP, reduction = "tsne",label=T)

plot11<-DimPlot(Ctrl, reduction = "umap",label=T)
plot12<-DimPlot(CP, reduction = "umap",label=T)
plot11+plot12

Annotation<-function(Merged_xxx){
meta<-Merged_xxx@meta.data
#turning Seurat matrix into SingleR matrix
SingleR_xxx<-GetAssayData(Merged_xxx,slot="data")
#general match to lablel dataset
Merged_xxx.match<-SingleR(test=SingleR_xxx,ref=MouseRNAseq,labels=MouseRNAseq$label.main)
table(Merged_xxx.match$labels,meta$seurat_clusters)
Merged_xxx@meta.data$labels <-Merged_xxx.match$labels
DimPlot(Merged_xxx,group.by=c("labels"),reduction = "umap")
ggplot(data=as.data.frame(table(Merged_xxx.match$labels,meta$seurat_clusters)),aes(x=Var1,y=Var2,size=Freq))+
geom_point(alpha=0.8)+scale_size(range=c(1,20))+theme_bw()+labs(x="SingleR",y="Seurat")
return(Merged_xxx)
}

#run the annotation function
Merged_Ctrl_Annotated<-Annotation(Ctrl)
Merged_CP_Annotated<-Annotation(CP)

#plot after annotation
'
cluster_color<-c("Adipocytes"="#655548","B cells"="#9F6763","Dendritic cells"="#F9CB19","Endothelial cells"="#EE5A9B","Epithelial cells"="#C392BC",
"Fibroblasts"="#8DC591","Granulocytes"="#C8A99B","Hepatocytes"="#BE9C6E","Macrophages"="#F58825","Microglia"="#EE5125","Monocytes"="#B61F8A",
"NK cells"="8770AC","T cells"="542E88")
'
cluster_color<-c("B cells"="#F9CB19","Dendritic cells"="#124D98","Endothelial cells"="#EE5A9B",
"Fibroblasts"="#F59084","Granulocytes"="#C2C287","Hepatocytes"="#BE9C6E","Macrophages"="#F58825","Microglia"="#EFEF9D","Monocytes"="#B61F8A",
"NK cells"="#8770AC","T cells"="#C7212A","Adipocytes"="#655548","Epithelial cells"="#D8BC97")

plot13<-DimPlot(Merged_Ctrl_Annotated, reduction = "umap",group.by=c("labels"),pt.size=0.5,label=T)+scale_color_manual(values=cluster_color)
plot14<-DimPlot(Merged_CP_Annotated, reduction = "umap",group.by=c("labels"),pt.size=1,label=T)+scale_color_manual(values=cluster_color)
plot13+plot14

save(Merged_Ctrl_Annotated,file="RData/Merged_Ctrl_Annotated.RData")
save(Merged_CP_Annotated,file="RData/Merged_CP_Annotated.RData")


###############################plot the marker genes#################################

#find the top 30 variable features and filter 5 from each cluster
top30_Ctrl <- head(VariableFeatures(Merged_Ctrl_Annotated), 30)
#show the head of the variable features
head(top30_Ctrl)
Merged_Ctrl_Annotated.markers <- FindAllMarkers(Merged_Ctrl_Annotated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Idents(Merged_Ctrl_Annotated)<-'labels'
Merged_Ctrl_Annotated.markers %>% group_by(gene) %>% top_n(n = 5, wt = avg_log2FC) %>% ungroup() %>% dplyr::filter(gene %in% top30_Ctrl)

top30_CP <- head(VariableFeatures(Merged_CP_Annotated), 30)
#show the head of the variable features
head(top30_CP)
Merged_CP_Annotated.markers <- FindAllMarkers(Merged_CP_Annotated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Idents(Merged_CP_Annotated)<-'labels'
Merged_CP_Annotated.markers %>% group_by(gene) %>% top_n(n = 5, wt = avg_log2FC) %>% ungroup() %>% dplyr::filter(gene %in% top30_CP)

##################
save(Merged_Ctrl_Annotated.markers,file="RData/Merged_Ctrl_Annotated.markers.RData")
save(Merged_CP_Annotated.markers,file="RData/Merged_CP_Annotated.markers.RData")

load("RData/Merged_Ctrl_Annotated.RData")
load("RData/Merged_CP_Annotated.RData")
load("RData/Merged_Ctrl_Annotated.markers.RData")
load("RData/Merged_CP_Annotated.markers.RData")
################

library(ggsci)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(svglite)
library(Seurat) #version 4.4.0
library(dplyr) #helps to rename the column
library(Matrix) #Matrix 1.6.1 or there will be an as_cholmod_sparse error
#install.packages("Matrix_1.6-1.tar.gz",repos = NULL,type="source")
library(patchwork)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

#select the top 5 genes from each cluster
top5_Ctrl.markers <- Merged_Ctrl_Annotated.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)


top5_CP.markers <- Merged_CP_Annotated.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

#general heatmap
GeneralHeatmap<-function(pbmc,top5pbmc.markers){
DoHeatmap(pbmc,
features=top5pbmc.markers$gene,group.colors=unname(cluster_color))+
scale_colour_npg()+scale_fill_gradient2(low='#0099CC',mid='white',high='#CC0033',nam='Z-score')
}

#cluster heatmap
ClusterHeatmap<-function(pbmc,top5pbmc.markers){

all_cell_exp<-FetchData(pbmc,vars=top5pbmc.markers$gene,slot='data')%>% t()
metadata<-pbmc@meta.data%>%arrange(labels)#help to order the cell expression, so that the heatmap is in ladder order
#metadata$cluster_name<-paste('cluster',metadata$seurat_clusters,sep='_')
#order cell
all_cell_exp<-all_cell_exp[,rownames(metadata)]
#Z-score
all_cell_exp<-t(scale(t(all_cell_exp),scale = T,center = T))
#check
allht<-t(scale(t(all_cell_exp),scale = T,center = T))
#color
col_fun=colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))
#top annotation
table(metadata$labels)
anno_col1<-brewer.pal(8,"Paired")
names(anno_col1)<-paste(metadata$labels,1:10,sep='_')
#names(anno_col)<-paste('cluster',1:10,sep='_')
#make annotation
column_ha<-HeatmapAnnotation(cluster=metadata$labels,
                            col = list(labels = anno_col1),border=T)
#plot
Heatmap(allht,
        name = "Z-score",
        cluster_columns = F,cluster_rows = F,
        row_title = "Cluster top 5 Marker genes",
        column_title = "Cell Type",
        row_names_gp = gpar(fontface = 'italic',fontsize = 10),
        show_column_names=F,
        row_names_side = 'left',
        border = T,
        rect_gp = gpar(col = "white", lwd = 1),
        #make the cluster color
        column_names_side = 'top',
        column_names_rot = 45,
        column_split=metadata$cluster_name,
        top_annotation = column_ha,
        col = col_fun
        #ggsave("MarkerHeat.svg",width = 10,height = 10)
        #save a pdf for the plot
        )
}
#Marker heatmap
MarkerHeatmap<-function(pbmc,top5pbmc.markers,filename,save_svg=FALSE,save_png=FALSE,save_pdf=FALSE){
cluster_color <- c("RNA.B.cells"="#655548","RNA.Dendritic.cells"="#F9CB19","RNA.Endothelial.cells"="#EE5A9B","RNA.Fibroblasts"="#F59084","RNA.Granulocytes"="#C2C287",
"RNA.Hepatocytes"="#BE9C6E","RNA.Macrophages"="#228B22","RNA.Microglia"="#EFEF9D","RNA.Monocytes"="#4682B4","RNA.NK.cells"="#124D98",
 "RNA.T.cells"="#C7212A", "RNA.Adipocytes"="#542E88","RNA.Epithelial.cells"="#8770AC")
# get top 5 genes
mean_gene_exp <- AverageExpression(pbmc,
                                   features = top5pbmc.markers$gene,
                                   group.by = 'labels',
                                   slot = 'data') %>%
                                    data.frame() %>%
                                    as.matrix()

# add colnames
scaled_expression<-t(scale(t(mean_gene_exp),scale = T,center = T))
gene_color=colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))
cluster_number<-(length(colnames(scaled_expression)))

annotation_color <- cluster_color[colnames(scaled_expression)]
names(annotation_color)<-colnames(scaled_expression)
heatmap_setting=HeatmapAnnotation(cluster=colnames(scaled_expression),
                            col = list(cluster = annotation_color),border=T)#add color to labels 

# plot
ht<-Heatmap(scaled_expression,
        name = "Z-score",
        cluster_columns = F,cluster_rows = F,
        row_title = "Cluster top 5 Marker genes",
        column_title = "Cell Type",
        row_names_gp = gpar(fontface = 'italic',fontsize = 10),
        row_names_side = 'left',
        border = T,
        rect_gp = gpar(col = "white", lwd = 1),
        #make the cluster color
        column_names_side = 'top',
        column_names_rot = 45,
        top_annotation = heatmap_setting,
        col = gene_color
        #ggsave("MarkerHeat.svg",width = 10,height = 10)
        #save a pdf for the plot
        )
if(save_svg){
  svg(filename,width=10,height=20) #only turn on when plot the general heatmap
  draw(ht)#only turn on when plot the general heatmap
  dev.off() ##only turn on when plot the general heatmap
}
else if(save_png){
  png(filename,width=10,height=20) #only turn on when plot the general heatmap
  draw(ht)#only turn on when plot the general heatmap
  dev.off() ##only turn on when plot the general heatmap
}
else if(save_pdf){
  pdf(filename,width=10,height=20) #only turn on when plot the general heatmap
  draw(ht)#only turn on when plot the general heatmap
  dev.off() ##only turn on when plot the general heatmap
}
else{
  draw(ht)
}
}

plot15<-GeneralHeatmap(Merged_Ctrl_Annotated,top5_Ctrl.markers)
plot16<-GeneralHeatmap(Merged_CP_Annotated,top5_CP.markers)
plot15+plot16

plot17<-MarkerHeatmap(Merged_Ctrl_Annotated,top5_Ctrl.markers)
plot18<-MarkerHeatmap(Merged_CP_Annotated,top5_CP.markers)
plot17
plot18


###################restart from the Merged ####################
#if the UMAP and heatmap looks fine, then merge the two datasets

Merged_Ctrl_Annotated$group<-"Ctrl"
Merged_CP_Annotated$group<-"CP"
Merged<-merge(Merged_Ctrl_Annotated,Merged_CP_Annotated)
#preprocess the merged data
Merged_processed<-Pre_process(Merged)
Merged_debatched<-debatch(Merged_processed)
top30_Merged <- head(VariableFeatures(Merged_debatched), 30)
head(top30_Merged)
Merged_debatched.markers <- FindAllMarkers(Merged_debatched, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Idents(Merged_debatched)<-'labels'
Merged_debatched.markers %>% group_by(gene) %>% top_n(n = 5, wt = avg_log2FC) %>% ungroup() %>% dplyr::filter(gene %in% top30_Merged)

top5_Merged_debatched.markers <- Merged_debatched.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
#annotate the merged data


Merged_final<-Annotation(Merged_debatched)
DimPlot(Merged_final,group.by=c("labels"),reduction="umap",pt.size=0.5,label=T)+scale_color_manual(values=cluster_color)
Merged_final.markers<-Merged_debatched.markers
top5_Merged_final.markers<-top5_Merged_debatched.markers

save(Merged_final,file="RData/Merged_final.RData")
save(Merged_final.markers,file="RData/Merged_final.markers.RData")
save(top5_Merged_final.markers,file="RData/top5_Merged_final.markers.RData")
###############decoration in UMAP #################################
load("RData/Merged_final.RData")
load("RData/Merged_final.markers.RData")

#predine the dim and colorz
dims<-8

#plot the general Markerheatmap
top3_Merged_final.markers <- Merged_final.markers %>%x
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC)

Final_Markerheatmap<-MarkerHeatmap(Merged_final,top3_Merged_final.markers,"Final_Markerheatmap.pdf")#,save_pdf=T
#save the plot in svg

cluster_color<-c("B cells"="#F9CB19","Dendritic cells"="#124D98","Endothelial cells"="#EE5A9B",
"Fibroblasts"="#F59084","Granulocytes"="#C2C287","Hepatocytes"="#BE9C6E","Macrophages"="#F58825","Microglia"="#EFEF9D","Monocytes"="#B61F8A", 
"NK cells"="#8770AC","T cells"="#C7212A","Adipocytes"="#655548","Epithelial cells"="#D8BC97")
cluster_color2<-c("B cells"="#F9CB19","Dendritic cells"="#124D98","Endothelial cells"="#EE5A9B",
"Fibroblasts"="#F59084","Granulocytes"="#C2C287","Hepatocytes"="#BE9C6E","Macrophages"="#563F2E","Microglia"="#78C2C4","Monocytes"="#C73E3A", 
"NK cells"="#A35E47","T cells"="#B47157","Adipocytes"="#2B5F75","Epithelial cells"="#D8BC97")
cluster_color3 <- c("B cells"="#655548","Dendritic cells"="#F9CB19","Endothelial cells"="#EE5A9B","Fibroblasts"="#F59084","Granulocytes"="#C2C287",
"Hepatocytes"="#BE9C6E","Macrophages"="#228B22","Microglia"="#EFEF9D","Monocytes"="#4682B4","NK cells"="#124D98",
 "T cells"="#C7212A", "Adipocytes"="#542E88","#8770AC")
library(ggunchull)
library(tidydr)
library(ggrepel)
#plot the general UMAP

#Final_UMAP<-DimPlot(Merged_final,split.by='group',reduction="umap",cols=cluster_color,pt.size=2)
Final_UMAP<-DimPlot(Merged_final,reduction="umap",cols=cluster_color,pt.size=2)
Merged_final_data<-Merged_final@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=Merged_final@meta.data$labels)
axis<-Merged_final_data %>%
  group_by(celltype) %>%
   dplyr::summarise(
      UMAP_1=median(UMAP_1),
      UMAP_2=median(UMAP_2))

UMAP<-Final_UMAP+geom_label_repel(data=axis,aes(x=UMAP_1,y=UMAP_2,label=celltype),fontface="bold",point.padding=unit(0.3,"lines"))+
theme_dr(xlength=0.3,ylength=0.3,arrow=grid::arrow(length=unit(0.2,"inches"),type="closed"))+
theme(panel.grid=element_blank())
#save the UMAP
ggsave("general_UMAP.pdf",plot=UMAP,width=10,height=10)
ggsave("general_UMAP.svg",plot=UMAP,width=10,height=10)

###########################plot the selected celltype#####################################
selected_celltype<-c("B cells","Dendritic cells","Macrophages","Monocytes","NK cells","T cells")

subset_final<-subset(Merged_final,subset=labels %in% selected_celltype)
Selected_UMAP_VS<-DimPlot(subset_final,reduction="umap",cols=cluster_color3,pt.size=1,group.by="labels",split.by="group")

#plot the  anntotation labels 
subset_final_data<-subset(Merged_final_data,celltype %in% selected_celltype)
axis<-subset_final_data %>%
  group_by(celltype) %>%
   dplyr::summarise(
      UMAP_1=median(UMAP_1),
      UMAP_2=median(UMAP_2)
      )
selected_UMAP<-Selected_UMAP_VS+geom_label_repel(data=axis,aes(x=UMAP_1,y=UMAP_2,label=celltype),fontface="bold",point.padding=unit(0.3,"lines"))+scale_fill_manual(values=cluster_color3)+
theme(legend.text=element_text(color="white"))+theme_dr(xlength=0.3,ylength=0.3,arrow=grid::arrow(length=unit(0.2,"inches"),type="closed"))+theme(panel.grid=element_blank())
ggsave("selected_UMAP.pdf",plot=selected_UMAP,width=16,height=10)+scale_color_manual(values=cluster_color3)
#draw the celltype ratio barplot
selected_ratio_bar<-ggplot(subset_final_data,aes(x=celltype,fill=celltype))+
  geom_bar()+
  facet_wrap(~group,scales="free")+
  scale_fill_manual(values=cluster_color3)+
  theme(legend.position="none",axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(x="Cell type",y="Cell number")+
  theme_dr(xlength=0.3,ylength=0.3,arrow=grid::arrow(length=unit(0.2,"inches"),type="closed"))+theme(panel.grid=element_blank())+
  coord_flip()+
  geom_text(aes(label=after_stat(count)),stat="count",vjust=-0.25)+
  theme(plot.margin=margin(1,1,1,7,"cm"))
ggsave("UMAP with counts.pdf",plot=(selected_UMAP+selected_ratio_bar),width=30,height=10)

#+stat_ellipse(aes(x=UMAP_1,y=UMAP_2,fill=cluster_color),geom="polygon",linetype=2,alpha=0.05,show.legend=FALSE,level=0.93)+
#theme_dr(xlength=0.22,ylength=0.22,arrow=grid::arrow(length=unit(0.15,"inches"),type="closed"))+
#theme(panel.grid=element_blank())+stat_unchull(fill="white",alpha=0,show.legend = F,nsm=30,nbin=150,sfac=1.2)
#+scale_color_manual(values=cluster_color)

#add the circle outline(don`t delete)
'
UMAP=ggplot(Merged_final_data, aes(UMAP_1, UMAP_2, color=celltype,fill = celltype))+
  geom_point(size = 1) + #细胞亚群点的大小
  scale_color_manual(values=cluster_color)+
  theme(panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        #panel.grid = element_blank(),
        # panel.grid.major = element_blank(), #网格线
        # panel.grid.minor = element_blank(), #网格线
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"),#背景色
        legend.title = element_blank(), #去掉legend.title 
        legend.key=element_rect(fill='white'), 
        legend.text = element_text(size=14), #设置legend标签的大小
        legend.key.size=unit(0.6,'cm') )+
  stat_unchull(
    fill = "white",
    alpha = 0,
    show.legend = FALSE,
    nsm = 30,
    nbin = 150,
    sfac = 1.2)+
guides(fill= guide_legend(override.aes = list(size = 4)))+ #图例点的大小
  scale_color_manual(values=cluster_color)+
  stat_ellipse(aes(x = UMAP_1, y = UMAP_2, fill = celltype),
              geom = "polygon",#多边形
              linetype=2, # 椭圆，虚线
              alpha = 0.05,  #椭圆不透明度
              show.legend = FALSE, #去掉椭圆对应的图例
              level = 0.93)+ #置信区间
  guides(fill= guide_legend(override.aes = list(size = 4)))+ #图例点的大小
  scale_color_manual(values=cluster_color)
'
ggsave("UMAP.pdf",plot=UMAP,width=25,height=15)

NK_cells<-subset(subset_final,labels=="NK cells")
NK_cells_subcluster<-Annotation(NK_cells)
NK_cells_subcluster_UMAP<-DimPlot(NK_cells_subcluster, reduction = "umap",group.by=c("labels"),pt.size=0.5,label=T)+scale_color_manual(values=cluster_color3)
#plot the NK MarkerHeatmap
#features<-c("Gzma","Klrb1c","Ncr1","Klra8","Vps37b","Crem","Abcb1b","Bcl2l11","St6galnac3","Ms4a4b","Ccl5","Fcgr3")#Cell jounrnal CP
#features<-c("Klrc2","Mki67","Cx3cr1","Il32","Znf90","Nr4a3","Dnajb1","Nfkbia","Gnly","Xcl1")#NC journal
#features<-c("S1pr5","Prf1","Itga","Ctla2a","Gzmc","Xcl1","Klrg1","S1pr5","Cma1","Nkg7")#CellMarker 2.0
features_NK<-c("Ncr1","Dnajb1","Crem","Bcl2l11")#Validated NK markers omit "Abcb1b","Nr4a3"
#NK_Ctrl<-FeaturePlot(Merged_Ctrl_Annotated,features=features,min.cutoff=1,max.cutoff=3)#validated
#NK_CP<-FeaturePlot(Merged_CP_Annotated,features=features,min.cutoff=1,max.cutoff=3)#validated
NK_markers<-FeaturePlot(NK_cells,features=features_NK,min.cutoff=1,max.cutoff=3,split.by="group")
#save the plot in pdf
ggsave("NK_markers.png",plot=NK_markers,width=15,height=20)

#plot the T MarkerHeatmap
#features<-c("Cd3g","Icos","Il23r","Camk4","Zbtb16")#Cell jounrnal Ctrl
#features<-c("Nebl","Inpp4b","Il7r","St6galnac3","Cd7","Cd3g","Icos","Rora","Rgcc","Cntn1","Lingo2","Dscam","Large1","Cd3d","Lef1","Grap2")#NC journal
features_T_cells<-c("Icos","Zbtb16","Cntn1","Cd3g")
#T_Ctrl<-FeaturePlot(Merged_Ctrl_Annotated,features=features,min.cutoff=1,max.cutoff=3)#validated
#T_CP<-FeaturePlot(Merged_CP_Annotated,features=features,min.cutoff=1,max.cutoff=3)#validated
T_markers<-FeaturePlot(subset_final,features=features_T_cells,min.cutoff=1,max.cutoff=3,split.by="group")
T_ridge<-RidgePlot(subset_final,features=features_T_cells,ncol=2)
#ggsave("T_markers.tif",plot=T_markers,width=15,height=20,device="tiff")
ggsave("T_markers.png",plot=T_markers,width=15,height=20)

Monocytes<-subset(subset_final,labels=="Monocytes")
Monocytes_subcluster<-Annotation(Monocytes)
Monocytes_subcluster_UMAP<-DimPlot(Monocytes_subcluster, reduction = "umap",group.by=c("labels"),pt.size=0.5,label=T)+scale_color_manual(values=cluster_color3)
features_Monocytes_Macrophages<-c("Fcgr2b","Cd68")#From CellMarker 2.0,complied with SinglR
features_Monocytes_Macrophages<-c("Chil3")#from Cell journal,complied with SinglR
#T_Ctrl<-FeaturePlot(Merged_Ctrl_Annotated,features=features,min.cutoff=1,max.cutoff=3)#validated
#T_CP<-FeaturePlot(Merged_CP_Annotated,features=features,min.cutoff=1,max.cutoff=3)#validated
Monocytes_markers<-FeaturePlot(subset_final,features=features_Monocytes_Macrophages,min.cutoff=1,max.cutoff=3,split.by="group",pt.size=0.8)
#ggsave("T_markers.tif",plot=T_markers,width=15,height=20,device="tiff")
ggsave("Monocytes_markers.png",plot=Monocytes_markers,width=24,height=12)

features_DC<-c("Itgax","Clecl9a","Batf3","Bst2")
DC_markers<-FeaturePlot(subset_final,features=features_DC,min.cutoff=1,max.cutoff=3,split.by="group",pt.size=0.8)
#ggsave("T_markers.tif",plot=T_markers,width=15,height=20,device="tiff")
ggsave("DC_markers.png",plot=Monocytes_markers,width=24,height=12)

#plot the MarkerHeatmap of the selected celltype
subset_final.markers <- FindAllMarkers(subset_final, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#plot the general Markerheatmap
top5_subset_final.markers <- subset_final.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

Final_selected_Markerheatmap<-MarkerHeatmap(subset_final,top5_subset_final.markers,"Final_selected_Markerheatmap.svg",,save_svg=T)

#######################cluster in subcluster of selected celltype####################
Macrophages<-subset(subset_final,labels=="Macrophages")
Macrophage_subcluster<-Annotation(Macrophages)
Macrophages_subcluster_UMAP<-DimPlot(Macrophage_subcluster, reduction = "umap",group.by=c("labels"),pt.size=0.5,label=T)+scale_color_manual(values=cluster_color3)
cellpred<-SingleR(test=GetAssayData(Macrophage_subcluster,slot="data"),ref=MouseRNAseq,labels=MouseRNAseq$label.main)
celltype=data.frame(ClusterID=rownames(cellpred),celltype=cellpred$labels,stringsFactors=F)
write.csv(celltype,"Macrophages_subcluster.csv")
for (i in 1:nrow(celltype)){
  subset_final@meta.data[which(subset_final@meta.data$seurat_clusters==celltype$ClusterID[i],"celltype")]=celltype$celltype[i]

}


