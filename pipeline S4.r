
#now you have the barcodes.tsv, features.tsv, matrix.mtx after running the cellranger

######remember using "str" often to check the structure whenever come across a new data #######

####################################build/activate the virtual environment###################

#conda create -n bio r-base
conda activate bio
#enter R
R

install.packages("Matrix_1.6-1.tar.gz",repos=NULL,type="source")#install locally
library(Matrix)
install.packages("SeuratObject_4.1.4.tar.gz",repos=NULL,type="source")
install.packages("Seurat")# This would help to install most required dependent packages even when Seurat-V5 will not be installed successfully
install.packages("remotes")#remotes was the simplified "devtools",help to install packages that are under developing(not on CRAN)
remotes::install_github("r-lib/devtools")
remotes::install_github('igraph/rigraph')
remotes::install_github('TomKellyGenetics/leiden')
install.packages("Seurat_4.4.0.tar.gz",repos=NULL,type="source")
library(Seurat) #version 4.4.0
install.packages("harmony")
library(harmony)
install.packages("ggplot2_3.5.0.tar.gz",repos=NULL,type="source")
library(ggplot2)
install.packages("tidygraph")
library(igraph)
BiocManager::install("SingleR")
library(SingleR)
###########################################load the data##########################################
#import the the extracted datas and merge
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
dim(group)
###############filter#############
load("./RData/group.RData")
Pre_process<-function(data){
data<-CreateSeuratObject(counts=data,project='Mouse',min.cells=3,min.features=200)
data[["percent_mt"]]<-PercentageFeatureSet(data,pattern="^mt-")
#refine the data
data<-subset(data,subset=nFeature_RNA>200&nFeature_RNA<8000&percent_mt<5)#####change from[200,2000] Dec.16 to [200,8000] Mar.27
#nFeature_RNA means unique gene, dead or damaged cells shows low nFeatrue_RNA, doublets or polluted cells shows high nFeature_RNA
# VlnPlot(data,features=c("nFeature_RNA","nCount_RNA","percent_mt"),ncol=3)
#pca
data<-NormalizeData(data)
data<-FindVariableFeatures(data,selection.method="vst",nfeatures=2000)
all.genes<-rownames(data)
data<-ScaleData(data,features=all.genes)
data<-RunPCA(data,features=VariableFeatures(object=data))
ElbowPlot(data)#dims were decided from here, where the standard deviation came to the minimum,in our project,dims is 20
dims=20
FeatureScatter(data,feature1="nCount_RNA",feature2="nFeature_RNA")
FeatureScatter(data,feature1="nCount_RNA",feature2="percent_mt")
VizDimLoadings(data,reduction="pca",dims=1:dims)   
DimHeatmap(data,dims=1:dims,cells=800,balanced=TRUE)
data<-FindNeighbors(data,dims=1:dims)
data<-FindClusters(data,resolution=1.5,method="snn",algorithm=1)
#data<-RunTSNE(data,dims=1:dims)
data<-RunUMAP(data,dims=1:dims)
#plot the UMAP
#DimPlot(data,reduction="tsne",label=T,group.by="seurat_clusters")
DimPlot(data,reduction="umap",label=T,split.by="seurat_clusters")
#sometimes cannont run UMAP because as_cholmod_sparse not provided by package "Matrix"
return(data)
}
processed_group<-lapply(group,FUN=Pre_process) # processe the  data and save them in a list
#names(processed_group) <- paste0("processed_",names(processed_group)) #add the prefix "processed_" to the names of the list
save(processed_group,file="RData/processed_group_nFea8k.RData")
save(processed_group,file="RData/processed_group_dims25.RData")
save(processed_group,file="RData/processed_group_dims20_nFea8k.RData")
##########################de doublet cells###############################
load("RData/processed_group_dims.RData")
load("RData/processed_group_dims25.RData")
library(DoubletFinder)#version 2.0.3 
#paramSweep_v3.R and DoubletFinder_v3.R have been modified ,where meta@data$RNA@counts is changed to GetAssayData(meta, "counts")
#remotes::install_github("chris-mcginnis-ucsf/DoubletFinder") remotely
#then remove.packages("DoubletFinder"),this shall only leave dependencies,and finally
#install.packages("./DoubletFinder-master",repos=NULL,type="source") locally
#make sure DoubletFinder has been installed in defaulted path
library(patchwork)
Find_doublet <- function(processed_xxx) {
sweep.res.list <-paramSweep_v3(processed_xxx, PCs=1:dims, sct = FALSE, num.cores=1)#seu@assays$RNA@counts to Seurat::GetAssayData(seu, "counts")
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() #choose the best pK
DoubletRate = ncol(processed_xxx) * 8 * 1e-6 #0.027296
nExp_poi <- round(DoubletRate * ncol(processed_xxx)) #93
#homotypic.prop <- modelHomotypic(data, PCs = 1:dims, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
processed_xxx <-doubletFinder_v3(processed_xxx, PCs = 1:dims, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
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
save(de2_group,file="./RData/de2_group_dims20_nFea8k.RData")

#################################################Merge###########################
#merge the three de2_group
load("RData/de2_group.RData")
load("RData/de2_group_dims20_nFea8k.RData")
#name the subgroup and regroup as CP and Ctrl
de2_group[[1]]@meta.data$orig.ident <- "CP1"
de2_group[[2]]@meta.data$orig.ident <- "CP2"
de2_group[[3]]@meta.data$orig.ident <- "CP3"
de2_group[[4]]@meta.data$orig.ident <- "Ctrl1"
de2_group[[5]]@meta.data$orig.ident <- "Ctrl2"
de2_group[[6]]@meta.data$orig.ident <- "Ctrl3"

Merged_CP<- merge(de2_group[[1]], y = list(de2_group[[2]], de2_group[[3]]))#sum all the cells [0327]
Merged_Ctrl<- merge(de2_group[[4]], y = list(de2_group[[5]], de2_group[[6]]))

save(Merged_CP,file="RData/Merged_CP_n800_rs0.8.RData")
save(Merged_Ctrl,file="RData/Merged_Ctrl_n800_rs0.8.RData")
##############################de batch effect####################################
load("RData/Merged_CP.RData")
load("RData/Merged_Ctrl.RData")

library(harmony) #version 1.2.0
#install.packages("harmony")
dims=20
debatch<-function(Merged_xxx){
Merged_xxx<-NormalizeData(Merged_xxx)
#Merged_xxx<-FindVariableFeatures(Merged_xxx)
all.genes_xxx<-rownames(Merged_xxx)
Merged_xxx<-ScaleData(Merged_xxx,features=all.genes_xxx)
Merged_xxx<-RunPCA(Merged_xxx,npcs=50,verbose=FALSE)
Merged_xxx<-RunHarmony(Merged_xxx,"orig.ident")
#Merged_xxx<-RunTSNE(Merged_xxx,dims=1:dims,reduction="harmony")
Merged_xxx<-RunUMAP(Merged_xxx,dims=1:dims,reduction="harmony")#error because of Matrix package(as_cholmod_sparse)
Merged_xxx<-FindNeighbors(Merged_xxx,dims=1:dims,reduction="harmony")
Merged_xxx<-FindClusters(Merged_xxx,resolution=1.8,reduction="harmony")
DimPlot(Merged_xxx,reduction="umap",label=T,group.by="orig.ident")
#DimPlot(Merged_xxx,reduction="tsne",label=T,group.by="orig.ident")
return(Merged_xxx)
}
Merged_Ctrl_debatched<-debatch(Merged_Ctrl)
Merged_CP_debatched<-debatch(Merged_CP)

Cluster_tsne_Ctrl<-DimPlot(Merged_Ctrl_debatched,reduction='umap',label=T,group.by='seurat_clusters')
Cluster_tsne_CP<-DimPlot(Merged_CP_debatched,reduction='umap',label=T,group.by='seurat_clusters')
Cluster_tsne_Ctrl+Cluster_tsne_CP

save(Merged_CP_debatched,file="RData/Merged_CP_debatched.RData")
save(Merged_Ctrl_debatched,file="RData/Merged_Ctrl_debatched.RData")
save(Merged_CP_debatched,file="RData/Merged_CP_debatched_nFea8k.RData")
save(Merged_Ctrl_debatched,file="RData/Merged_Ctrl_debatched_nFea8k.RData")
save(Merged_CP_debatched,file="RData/Merged_CP_debatched_dims36_res1.8.RData")#[0316]       [0327]#change filename from ....res1.5 to ...res1.8,align with Pre_process and debatch
save(Merged_Ctrl_debatched,file="RData/Merged_Ctrl_debatched_dims36_res1.8.RData")#[0316]

############merge CP and Ctrl and pre_process ##############################
Pre_process2<-function(data){
data<-NormalizeData(data)
data<-FindVariableFeatures(data,selection.method="vst",nfeatures=2000)
all.genes<-rownames(data)
data<-ScaleData(data,features=all.genes)
data<-RunPCA(data,features=VariableFeatures(object=data))
ElbowPlot(data)#dims were decided from here, where the standard deviation came to the minimum,in our project,dims is 20
dims=20
FeatureScatter(data,feature1="nCount_RNA",feature2="nFeature_RNA")
FeatureScatter(data,feature1="nCount_RNA",feature2="percent_mt")
VizDimLoadings(data,reduction="pca",dims=1:dims)   
DimHeatmap(data,dims=1:dims,cells=800,balanced=TRUE)
data<-FindNeighbors(data,dims=1:dims)
data<-FindClusters(data,resolution=0.9,method="snn",algorithm=1)#clusters coudl reach at 45 if res stays at 1.5
#data<-RunTSNE(data,dims=1:dims)
data<-RunUMAP(data,dims=1:dims)
#plot the UMAP
#DimPlot(data,reduction="tsne",label=T,group.by="seurat_clusters")
DimPlot(data,reduction="umap",label=T,split.by="seurat_clusters")
#sometimes cannont run UMAP because as_cholmod_sparse not provided by package "Matrix"
return(data)
}
Merged_debatched<-Pre_process2(merge(Merged_Ctrl_debatched,Merged_CP_debatched))
Merged_debatched@meta.data$group<-gsub("Ctrl[1-3]","Ctrl",Merged_debatched@meta.data$orig.ident)
Merged_debatched@meta.data$group<-gsub("CP[1-3]","CP",Merged_debatched@meta.data$group)
Merged_debatched@meta.data$group<-factor(Merged_debatched@meta.data$group,levels=c("Ctrl","CP"))

save(Merged_debatched,file="RData/Merged_debatched.RData")
save(Merged_debatched,file="RData/Merged_debatched_dims36_res1.8.RData")#[0316]
save(Merged_debatched,file="RData/Merged_debatched_dims20_nFea8k.RData")#[0327]
############################### Automated Annotation with SingleR (only for reference)#################################
load("RData/Merged_Ctrl_debatched.RData")
load("RData/Merged_CP_debatched.RData")
load("RData/Merged_debatched.RData")
load("RData/Merged_debatched_dims20_nFea8k.RData")
load("RData/Merged_debatched_dims36_res1.8.RData")

library(SingleR) #version 2.4.1 depended on Matrix 1.6-1,igraph,tidygraph
#install.packages("SingleR_2.4.1.tar.gz",repos = NULL,type="source")
library(celldex) #version 1.10.1
library(pheatmap)
library(RColorBrewer)
library(ggsci)
library(ggplot2)#depend on package tidyverse
#import mouse gene marker library manually(cuz automated always failed)
MouseRNAseq<-readRDS("MouseRNAseqData.rds") #or
#ImmuneGene<-get(load("ref_Mouse_imm.RData"))#a much richer datasets
Annotation<-function(Merged_xxx){
meta<-Merged_xxx@meta.data
#turning Seurat matrix into SingleR matrix
SingleR_xxx<-GetAssayData(Merged_xxx,slot="data")
#general match to lablel dataset
Merged_xxx.match<-SingleR(test=SingleR_xxx,ref=MouseRNAseq,labels=MouseRNAseq$label.main)
#Merged_xxx.match<-SingleR(test=SingleR_xxx,ref=ImmuneGene,labels=ImmuneGene$label.fine,assay.type.test="logcounts",assay.type.ref="logcounts")
table(Merged_xxx.match$labels,meta$seurat_clusters)
Merged_xxx@meta.data$labels <-Merged_xxx.match$labels
DimPlot(Merged_xxx,group.by=c("labels"),reduction = "tsne")
ggplot(data=as.data.frame(table(Merged_xxx.match$labels,meta$seurat_clusters)),aes(x=Var1,y=Var2,size=Freq))+
geom_point(alpha=0.8)+scale_size(range=c(1,20))+theme_bw()+labs(x="SingleR",y="Seurat")
return(Merged_xxx)
}
#run the annotation function
Merged_Annotated<-Annotation(Merged_debatched)
'
cluster_color<-c("Adipocytes"="#655548","B cells"="#9F6763","Dendritic cells"="#F9CB19","Endothelial cells"="#EE5A9B","Epithelial cells"="#C392BC",
"Fibroblasts"="#8DC591","Granulocytes"="#C8A99B","Hepatocytes"="#BE9C6E","Macrophages"="#F58825","Microglia"="#EE5125","Monocytes"="#B61F8A",
"NK cells"="8770AC","T cells"="542E88")
'
cluster_color<-c("B cells"="#F9CB19","Dendritic cells"="#124D98","Endothelial cells"="#EE5A9B",
"Fibroblasts"="#F59084","Granulocytes"="#C2C287","Hepatocytes"="#BE9C6E","Macrophages"="#F58825","Microglia"="#EFEF9D","Monocytes"="#B61F8A",
"NK cells"="#8770AC","T cells"="#C7212A","Adipocytes"="#655548","Epithelial cells"="#D8BC97") #organge color dominated
cluster_color<-c("NK cells"="#836621","T cells"="#c9c941","Monocytes"="#ae4ef3","Dendritic cells"="#5ac8c9",
"B cells"="#7a73f4","Endothelial cells"="#5cc86e","Fibroblasts"="#e83123","Macrophages"="#F58825","Adipocytes"="#655548","Epithelial cells"="#D8BC97",
"Granulocytes"="#C2C287","Hepatocytes"="#BE9C6E","Microglia"="#EFEF9D")#brown and silver color dominated

#Merged_debatched_tsne<-DimPlot(Merged_debatched,reduction = "tsne",label=T)
Merged_debatched_umap<-DimPlot(Merged_debatched,reduction = "umap",label=T)
ggsave("Merged_debatched_umap.pdf",plot=Merged_debatched_umap,width=10,height=10)
#Merged_debatched_tsne+Merged_debatched_umap
#ggsave("Merged_debatched_tsne_VS_umap.pdf",plot=Merged_debatched_tsne+Merged_debatched_umap,width=20,height=10)

Merged_Annotated_tsne<-DimPlot(Merged_Annotated,reduction = "tsne",group.by=("labels"),pt.size=1,label=F)+
scale_color_manual(values=cluster_color)
Merged_Annotated_umap<-DimPlot(Merged_Annotated,reduction = "umap",group.by=("labels"),pt.size=1,label=F)+
scale_color_manual(values=cluster_color)
Merged_Annotated_tsne+Merged_Annotated_umap
ggsave("Merged_Annotated_tsne_VS_umap.pdf",plot=Merged_Annotated_tsne+Merged_Annotated_umap,width=20,height=10)

save(Merged_Annotated,file="RData/Merged_Annotated_dims25.RData")
save(Merged_Annotated,file="RData/Merged_Annotated.RData")

###############################plot the marker genes#################################
load("RData/Merged_Annotated.RData")

top30 <- head(VariableFeatures(Merged_debatched), 30)
#show the head of the variable features
head(top30)

Merged_debatched.markers <- FindAllMarkers(Merged_debatched, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Merged_debatched.markers %>% group_by(gene) %>% top_n(n = 5, wt = avg_log2FC) %>% ungroup() %>% dplyr::filter(gene %in% top30)

save(Merged_Annotated.markers,file="RData/Merged_Annotated.markers.RData")
save(Merged_Annotated.markers,file="RData/Merged_Annotated_dims25.markers.RData")

load("RData/Merged_Annotated_dims25.markers.RData")
load("RData/Merged_Annotated_dims25.RData")

#############plot all markers heatmap(help to define cluster as celltype)####################
install.packages(c("ggsci","circlize","RColorBrewer","ComplexHeatmap","patchwork","pheatmap","ggplot2"))
library(ggsci)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
library(patchwork)
library(pheatmap)
library(ggplot2)
#select the top 5 genes from each cluster
top5_Merged_debatched.markers <- Merged_debatched.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

top5_Merged_Annotated.markers <- Merged_Annotated.markers %>%
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
all_cell_exp<-FetchData(pbmc,vars=top5pbmc.markers$gene,slot='data')
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
anno_col<-brewer.pal(12,"Paired")
names(anno_col)<-unique(metadata$labels)
#names(anno_col)<-paste('cluster',1:10,sep='_')
#make annotation
column_ha<-HeatmapAnnotation(cluster=metadata$labels,
                            col = list(labels = anno_col),border=T)
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
cluster_color <- c("RNA.B.cells"="#7a7df4","RNA.Dendritic.cells"="#5ac8c9","RNA.Endothelial.cells"="#e83123","RNA.Fibroblasts"="#96c8f2","RNA.Granulocytes"="#C2C287",
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
top3_Merged_debatched.markers <- top5_Merged_debatched.markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC)
top3_Merged_debatched.markers<-top5_Merged_Annotated.markers
GeneralHeatmap<-GeneralHeatmap(Merged_debatched,top3_Merged_debatched.markers)
MarkerHeatmao<-MarkerHeatmap(Merged_Annotated,top5_Merged_Annotated.markers)
##save(Merged_final,file="RData/Merged_final.RData")
save(Merged_Annotated,file="RData/Merged_Annotated.RData")
save(Merged_Annotated.markers,file="RData/Merged_Annotated.markers.RData")
save(top5_Merged_Annotated.markers,file="RData/top5_Merged_Annotated.markers.RData")


###############decoration in UMAP #################################
load("RData/Merged_Annotated.RData")
load("RData/top5_Merged_Annotated.markers.RData")

library(ggunchull)
#remotes::install_github("sajuukLyu/ggunchull",type="source")
library(tidydr)
#install.packages("tidydr")
library(ggrepel)
#plot the general UMAP/TSNE

#predine the dim 
dims<-20

#plot the general Markerheatmap
top3_Merged_Annotated.markers <- top5_Merged_Annotated.markers %>%
  group_by(cluster) %>%
  top_n(n = 3, wt = avg_log2FC)

'
Merged_Annotated$labels[Merged_Annotated$labels=="Microglia"]<-"Macrophages"
Merged_Annotated$labels[["Microglia"]]<-"Macrophages"
Cluster_MarkerHeatmap<-ClusterHeatmap(Merged_debathced,top3_Merged_Annotated.markers)
Final_MarkerHeatmap<-MarkerHeatmap(Merged_Annotated,top3_Merged_Annotated.markers,"xxx.pdf",save_pdf=T)

cluster_color<-cluster_color_manually
final_TSNE<-DimPlot(Merged_Annotated,reduction="tsne",cols=cluster_color,pt.size=1,label=T,group.by="labels")
final_UMAP<-DimPlot(Merged_Annotated,reduction="umap",cols=cluster_color,pt.size=1,label=T,group.by="labels")
final_TSNE+final_UMAP
"#Final_UMAP<-DimPlot(Merged_final,reduction="umap",cols=cluster_color,pt.size=2)
Merged_Annotated_data<-Merged_Annotated@reductions$tsne@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=Merged_Annotated@meta.data$labels)
axis<-Merged_Annotated_data %>%
  group_by(celltype) %>%
   dplyr::summarise(
      tSNE_1=median(tSNE_1),
      tSNE_2=median(tSNE2))
 
TSNE<-Final_TSNE+geom_label_repel(data=axis,aes(x=tSNE_1,y=tSNE_2,label=celltype),fontface="bold",point.padding=unit(0.3,"lines"))+
theme_dr(xlength=0.3,ylength=0.3,arrow=grid::arrow(length=unit(0.2,"inches"),type="closed"))+theme(panel.grid=element_blank())
#UMAP<-Final_UMAP+geom_label_repel(data=axis,aes(x=UMAP_1,y=UMAP_2,label=celltype),fontface="bold",point.padding=unit(0.3,"lines"))+
#theme_dr(xlength=0.3,ylength=0.3,arrow=grid::arrow(length=unit(0.2,"inches"),type="closed"))+
#theme(panel.grid=element_blank())

#save the UMAP
#ggsave("general_UMAP.pdf",plot=UMAP,width=10,height=10)
ggsave("TSNE_VS_UMAP.pdf",plot=final_TSNE+final_UMAP,width=20,height=10)
'
###########################identify possible subtypes by markers checking #####################################
library(ggplot2)
library(dplyr)

#DimPlot(Merged_Annotated,reduction="tsne",cols=cluster_color,pt.size=1,label=T,split.by="group")
#final_TSNE<-DimPlot(Merged_Annotated,reduction="tsne",cols=cluster_color,pt.size=2,label=T,split.by="group")
#ggsave("tSNE-2.pdf",plot=final_TSNE_VS,width=20,height=10)

cluster_UMAP<-DimPlot(Merged_debatched,reduction="umap",pt.size=1,label=T)
ggsave("cluster_umap.pdf",plot=cluster_UMAP,width=10,height=10)

cell_counts<-table(test@meta.data$group,test@meta.data$labels)
cell_counts_df<-as.data.frame(cell_counts)
total_counts<-aggregate(Count~Group,cell_counts_df,sum)
names(cell_counts_df)<-c("Group","CellType","Count")

celltype_bar<-ggplot(cell_counts_df,aes(x=Group,y=Count,fill=CellType))+
geom_bar(stat="identity")+theme_minimal()+labs(x="Group",y="Count",fill="Cell Type")+
theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())+scale_fill_manual(values=cluster_color)+
geom_text(data=total_counts,aes(label=Count,y=Count),vjust=-0.5)

ggsave("Celltype_bar.pdf",plot=celltype_bar,width=4,height=6)
'selected_celltype<-c("B cells","Dendritic cells","Macrophages","Monocytes","NK cells","T cells")

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
selected_ratio_area<-ggplot(subset_final_data,aes(x=celltype,fill=celltype))+
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
ggsave("UMAP.pdf",plot=UMAP,width=25,height=15)
'
##########check the marker distribution by cluster_UMAP, and define the celltype manually##############

#panmarker checking
features_cd45<-c("Ptprc")#expression in all leukocytes
cd45_markers<-FeaturePlot(Merged_debatched,features=features_cd45,min.cutoff=1,max.cutoff=3)
features_cd11b<-c("Itgam")#expression in all leukocytes
cd11b_markers<-FeaturePlot(Merged_debatched,features=features_cd11b,min.cutoff=1,max.cutoff=3)
cd45_markers+cd11b_markers #cd45 looks fine, cd11b sucks

#B_cells<-subset(Merged_Annotated,labels=="B cells")
#B_cells_subcluster_umap<-DimPlot(B_cells, reduction = "umap",pt.size=2,label=T)+scale_color_manual(values=cluster_color)
features_B<-c("Cd79a","Cd19","Ms4a1")
B_markers<-FeaturePlot(Merged_debatched,features=features_B,min.cutoff=1,max.cutoff=3)#cluster 17

#Endothelia_cells<-subset(Merged_Annotated,labels=="Endothelia")
#Endothelia_subcluster_umap<-DimPlot(Endothelia, reduction = "umap",pt.size=2,label=T)+scale_color_manual(values=cluster_color)
features_Endothelia<-c("Pecam1","Cldn5","Prox1")
Endothelia_markers<-FeaturePlot(Merged_debatched,features=features_Endothelia,min.cutoff=1,max.cutoff=3)#cluster 20

#Epithelia_cells<-subset(Merged_Annotated,labels=="Epithelia")
#Epithelia_subcluster_umap<-DimPlot(Epithelia, reduction = "umap",pt.size=2,label=T)+scale_color_manual(values=cluster_color)
features_Epithelia<-c("Krt18","Epcam","Krt8")
Epithelia_markers<-FeaturePlot(Merged_debatched,features=features_Epithelia,min.cutoff=1,max.cutoff=3)#cluster 25

#Fibroblast_cells<-subset(Merged_Annotated,labels=="Fibroblast")
#Fibroblast_subcluster_umap<-DimPlot(Fibroblast, reduction = "umap",pt.size=2,label=T)+scale_color_manual(values=cluster_color)
features_Fibroblast<-c("Col1a1","Col1a2","Dcn","Ccl19","Cxcl12")#"Ccl19","Cxcl12" from cellmarker
Fibroblast_markers<-FeaturePlot(Merged_debatched,features=features_Fibroblast,min.cutoff=1,max.cutoff=3)#cluster 29

#DC_cells<-subset(Merged_Annotated,labels=="Dendritic cells")
#DC_subcluster_umap<-DimPlot(DC_cells, reduction = "umap",pt.size=2,label=T)+scale_color_manual(values=cluster_color)
features_cDC<-c("Batf3","Clec9a","Ccr7","Ccl22")
cDC_markers<-FeaturePlot(Merged_debatched,features=features_cDC,min.cutoff=1,max.cutoff=3,pt.size=0.8)#4,2,15,6,23,5,18,32,30
features_pDC<-c("Siglech","Ly6d")
pDC_markers<-FeaturePlot(Merged_debatched,features=features_pDC,min.cutoff=1,max.cutoff=3,pt.size=0.8)#cluster 8,3,16,33,26,28,31
cDC_markers+pDC_markers
#NK_cells<-subset(Merged_Annotated,labels=="NK cells")
#plot the MarkerHeatmap of the selected celltype
#NK_cells_subcluster_tsne<-DimPlot(NK_cells, reduction = "tsne",split.by="group",pt.size=2,label=T)+scale_color_manual(values=cluster_color)
#plot the NK MarkerHeatmap
#features<-c("Gzma","Klrb1c","Ncr1","Klra8","Vps37b","Crem","Abcb1b","Bcl2l11","St6galnac3","Ms4a4b","Ccl5","Fcgr3")#Cell jounrnal CP
#features<-c("Klrc2","Mki67","Cx3cr1","Il32","Znf90","Nr4a3","Dnajb1","Nfkbia","Gnly","Xcl1")#NC journal
#features<-c("S1pr5","Prf1","Itga","Ctla2a","Gzmc","Xcl1","Klrg1","S1pr5","Cma1","Nkg7")#CellMarker 2.0
features_NK<-c("Ncr1","Dnajb1","Crem","Bcl2l11","Eomes","Klrb1c")#Validated NK markers omit "Abcb1b","Nr4a3"
#NK_Ctrl<-FeaturePlot(Merged_Ctrl_Annotated,features=features,min.cutoff=1,ma<-x.cutoff=3)#validated
#NK_CP<-FeaturePlot(Merged_CP_Annotated,features=features,min.cutoff=1,max.cutoff=3)#validated
NK_markers<-FeaturePlot(Merged_debatched,features=features_NK,min.cutoff=1,max.cutoff=3)#cluster 13,11,12,21,19,22,7,1,0,10
#save the plot in pdf
ggsave("NK_markers.pdf",plot=NK_markers,width=15,height=20)

#T_cells<-subset(Merged_Annotated,labels=="T cells")
#T_subcluster_umap<-DimPlot(T_cells, reduction = "umap",pt.size=2,label=T)+scale_color_manual(values=cluster_color)
features_T<-c("Cd3e","Cd3d")
T_markers<-FeaturePlot(Merged_debatched,features=features_T,min.cutoff=1,max.cutoff=3,pt.size=0.8)#cluster 13,19

#Neutrophil_cells<-subset(Merged_Annotated,labels=="Neutrophil")
#Neutrophil_subcluster_umap<-DimPlot(Neutrophil_cells, reduction = "umap",pt.size=2,label=T)+scale_color_manual(values=cluster_color)
features_Neutrophil<-c("S100a9","S100a8","Csf3r")
Neutrophil_markers<-FeaturePlot(Merged_debatched,features=features_Neutrophil,min.cutoff=1,max.cutoff=3,pt.size=0.8)#cluster 27

#Mono_Mac<-subset(Merged_Annotated,labels==c("Manocytes/Macrophages"))
#Mono_Mac_umap<-DimPlot(Mono_Mac, reduction = "umap",group.by=c("labels"),pt.size=0.5,label=T)+scale_color_manual(values=cluster_color)
features_Mono_Mac<-c("Lyz2","C1qa","C1qb","Xcr1")
Mono_Mac_markers<-FeaturePlot(Merged_debatched,features=features_Mono_Mac,min.cutoff=1,max.cutoff=3,pt.size=0.8)#14,24
ggsave("Monocytes_markers.png",plot=Monocytes_markers,width=24,height=12)

#Lti_cells<-subset(Merged_Annotated,labels=="Lti")
#Lti_subcluster_umap<-DimPlot(Lti_cells, reduction = "umap",pt.size=2,label=T)+scale_color_manual(values=cluster_color)
#features_Lti<-c("Ccr6","Kit","Cd4")
#Lti_markers<-FeaturePlot(Merged_debatched,features=features_Lti,min.cutoff=1,max.cutoff=3,pt.size=0.8)

#features_Ilc2<-c("Il1rl1","Rora","Gata3","Icos")
#Ilc2_markers<-FeaturePlot(Merged_debatched,features=features_Ilc2,min.cutoff=1,max.cutoff=3,pt.size=0.8)

features_Mast<-c("Ms4a2","Cpa3")
Mast_markers<-FeaturePlot(Merged_debatched,fe atures=features_Mast,min.cutoff=1,max.cutoff=3,pt.size=0.8)

#plot the MarkerHeatmap of the selected celltype
Merged_debatched.markers <- FindAllMarkers(Merged_debatched, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#plot the general Markerheatmap
top5_Merged_debatched.markers <- Merged_debatched.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
Merged_debatched.markers<-Merged_debatched.markers 
Markerheatmap<-MarkerHeatmap(Merged_debatched,top5_Merged_Annotated.markers,"Markerheatmap.pdf",,save_pdf=T)
Markerheatmap_res1_5<-MarkerHeatmap(Merged_debatched,top5_Merged_debatched.markers,"Markerheatmap.pdf",,save_pdf=T)

#########Now you have varified the celltype by checking markers##############
#######################cluster by markers manually########################

#Merged_debatched<-debatch(Merged_debatched)
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
#show the marker genes by cluster

#plot UMAP by varified celltype
Merged_Manually_Annotated<-Merged_debatched
#define the cluster with celltype manually
new.cluster.id<-c(
  "NK cells","NK cells","cDC","pDC","cDC","cDC","cDC","NK cells",
  "pDC","Mono/Mac","NK cells","NK cells","NK cells","T cells","Mono/Mac",
  "cDC","pDC","B cells","cDC","T cells","Endothelia","NK cells","NK cells","cDC","Mono/Mac",
  "Epithelia","pDC","Neutrophil","pDC","Fibrolasts","cDC","pDC","cDC",
  "pDC","Mast"
)

names(new.cluster.id)<-levels(Merged_Manually_Annotated)
Merged_new<-RenameIdents(Merged_Manually_Annotated,new.cluster.id)
Merged_Manually_Annotated$celltype<-Idents(Merged_new)

#plot the marker by celltype
cluster_color_manually<-c("NK cells"="#806420","T cells"="#c4c43f","Mono/Mac"="#912dc4","cDC"="#5ac3c4","Neutrophil"="#c22029",
"B cells"="#7770ee","Endothelia"="#124b94","Fibroblasts"="#e85897","Epithelia"="#d3B893","Mast"="#5ac36b","pDC"="#017ee0")

Manually_UMAP<-DimPlot(Merged_Manually_Annotated, reduction = "umap",group.by="celltype",pt.size=0.3,label=F)+scale_color_manual(values=cluster_color_manually)+theme(panel.grid=element_blank())

Manually_UMAP<-Manually_UMAP+theme_dr(xlength=0.3,ylength=0.3,arrow=grid::arrow(length=unit(0.2,"inches"),type="closed"))+theme(panel.grid=element_blank())

ggsave("Manually_UMAP.pdf",plot=Manually_UMAP,width=10,height=10)
save(Merged_Manually_Annotated,file="RData/Merged_Manually_Annotated_nfea8k.RData")
#########show bubble plot#####
features_all_markers1<-c("Cd79a","Ms4a1","Pecam1","Cldn5","Cd3d","Ncr1","Batf3","Clec9a","Ccr7","Ccl22","Siglech","Ly6d","Bcl2l11","Eomes","Klrb1c","Cd3e","Ms4a2","Cpa3")
dotplot1<-DotPlot(object=Merged_Manually_Annotated,features=features_all_markers1,group.by='celltype',dot.scale=9,cols=c("#d2d2d2","#8e3e2f"))+
coord_flip()+RotatedAxis()
features_all_markers1<-c()
cluster_marker1<-dotplot1+theme(axis.text.x=element_text(color=sapply(levels(Merged_Manually_Annotated$celltype),get_color)))#+scale_y_discrete(position="right") #color the celltype text

features_all_markers2<-c("Gzma","Gzmb","Klrb1c","Foxp3","Tbet","Tcf","Cd3e","Sox4","S100a9","S100a8","Fcgr2b","Lyz2","C1qa","C1qb","Il1rl1","Gata3","Krt18","Epcam","Krt8")
dotplot2<-DotPlot(Merged_Manually_Annotated,features=features_all_markers2,group.by='celltype',dot.scale=9,cols=c("#d2d2d2","#8e3e2f"))+
coord_flip()+RotatedAxis()
cluster_marker2<-dotplot2+theme(axis.text.x=element_text(color=sapply(levels(Merged_Manually_Annotated$celltype),get_color)))#+scale_y_discrete(position="right") #color the celltype text
cluster_marker1+cluster_marker2

ggsave("cluster_marker.pdf",plot=cluster_marker1+cluster_marker2,width=12,height=6)

######plot the celltype ratio################
Ctrl_Manually_Annotated<-subset(Merged_Manually_Annotated,subset=group=="Ctrl")
CP_Manually_Annotated<-subset(Merged_Manually_Annotated,subset=group=="CP")

Ctrl_Manually_Annotated_data<-Ctrl_Manually_Annotated@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=Ctrl_Manually_Annotated$celltype)
axis<-Ctrl_Manually_Annotated_data %>%
  group_by(celltype) %>%
   dplyr::summarise(
      UMAP_1=median(UMAP_1),
      UMAP_2=median(UMAP_2))

CP_Manually_Annotated_data<-CP_Manually_Annotated@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype=CP_Manually_Annotated$celltype)
axis<-CP_Manually_Annotated_data %>%
  group_by(celltype) %>%
   dplyr::summarise(
      UMAP_1=median(UMAP_1),
      UMAP_2=median(UMAP_2))

celltype_ratio_Ctrl<-ggplot(Ctrl_Manually_Annotated_data,aes(x=celltype,fill=celltype))+geom_bar(position="stack")+
scale_fill_manual(values=cluster_color_manually)+theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1,vjust=1))+
labs(x="Cell type",y="Number of cells")+theme(panel.grid=element_blank(),panel.background=element_blank())+
geom_text(aes(label=after_stat(count)),stat="count",vjust=-0.25)+theme(plot.margin=margin(1,1,1,7,unit="cm"))

celltype_ratio_CP<-ggplot(CP_Manually_Annotated_data,aes(x=celltype,fill=celltype))+geom_bar(position="stack")+
scale_fill_manual(values=cluster_color_manually)+theme(legend.position="none",axis.text.x=element_text(angle=45,hjust=1,vjust=1))+
labs(x="Cell type",y="Number of cells")+theme(panel.grid=element_blank(),panel.background=element_blank())+
geom_text(aes(label=after_stat(count)),stat="count",vjust=-0.25)+theme(plot.margin=margin(1,1,1,7,unit="cm"))

celltype_ratio_Ctrl+celltype_ratio_CP

ggsave("celltype_ratio_VS.pdf",plot=celltype_ratio_Ctrl+celltype_ratio_CP,width=13,height=6)

################Let`s Go Deeper!!!####################
#################Subtype cluster mining#####################
#####Find the most different celltypes between CP and Ctrl, then plot the Markermap and Violin of that celltype, identify the 
##### specific marker for each cluster and name them each as (marker positive) subcelltype, then relate their expression diffrence to vaccine

#####NK cells#####
selected_NK_cells<-subset(Merged_Manually_Annotated,subset=celltype %in% c("NK cells"))
DimPlot(selected_NK_cells,reduction="umap",pt.size=1,label=T,split.by="group",group.by="seurat_clusters")
selected_NK_cells_0.5.markers<-FindAllMarkers(selected_NK_cells,only.pos=TRUE,min.pcty=0.5,logfc.threshold=0.5)
selected_NK_cells_0.5.markers%>%
            group_by(cluster)%>%
            dplyr::filter(avg_log2FC>1)%>%
            slice_head(n=10)%>%
            ungroup()->top100
NK_marker_0.5<-DoHeatmap(selected_NK_cells,features=top100$gene)

#[0]:"Cma1"
#[2,6,7,8,10,12,15,32]:"Nfkb1","Zbtb46","Vegfa"
#[5,10]:"Nebl"
#[26,32]:"Oasl1","Oasl2","Ifit3","Ifit1"
#[34]:"Cd3e","Cd3g" 
NK_violin<-VlnPlot(selected_NK_cells,features=c("Ifit1","Ifit3","Cd3e","Cd3d","Hspa1a","Gzma","Gzmb","Ccl4","Ifng"),
split.by="group",ncol=3,add.noise=FALSE,pt.size=0)
ggsave("NK_violin.pdf",plot=NK_violin,width=12,height=)

NK_Manually_Annotated<-selected_NK_cells
#define the cluster with celltype manually
new.cluster.id<-c(
  "NK cells(Gzma/b+)","NK cells(Gzma/b+)","NK cells(Gzma/b+)","NK cells(Gzma/b+)","NK cells(Gzma/b+)","NK cells(Gzma/b+)",
  "NK cells(Gzma/b+)","NK cells(Gzma/b+)","NK cells(Gzma/b+)","NK cells(Gzma/b+)","NK cells(Gzma/b+)","NK cells(Ifit+)",
  "NK cells(Ifit+)"
)
names(new.cluster.id)<-levels(NK_Manually_Annotated)
Merged_new<-RenameIdents(NK_Manually_Annotated,new.cluster.id)
NK_Manually_Annotated$celltype<-Idents(Merged_new)

#plot the marker by celltype
cluster_color_NK<-c("NK cells(Gzma/b+)"="#CE6D39","NK cells(Ifit+)"="#F00E42","NK cells(Cd3e/d+)"="#FFEEE4")
NK_Manually_UMAP<-DimPlot(NK_Manually_Annotated, reduction = "umap",group.by="celltype",pt.size=0.3,label=F)+theme(panel.grid=element_blank())+
theme_dr(xlength=0.3,ylength=0.3,arrow=grid::arrow(length=unit(0.2,"inches"),type="closed"))+theme(panel.grid=element_blank())+scale_color_manual(values=cluster_color_NK)
ggsave("NK_Manually_UMAP.pdf",plot=NK_Manually_UMAP,width=12,height=12)

####pDC############
selected_pDC<-subset(Merged_Manually_Annotated,subset=celltype %in% c("pDC"))
DimPlot(selected_pDC,reduction="umap",pt.size=1,label=T,split.by="group",group.by="seurat_clusters")
selected_pDC.markers<-FindAllMarkers(selected_pDC,only.pos=TRUE)
selected_pDC.markers%>%
            group_by(cluster)%>%
            dplyr::filter(avg_log2FC>1)%>%
            slice_head(n=10)%>%
            ungroup()->top50       
pDC_marker<-DoHeatmap(selected_pDC,features=top50$gene)

VlnPlot(selected_pDC,features=c())

####Neutrophil############
selected_Neutrophil<-subset(Merged_Manually_Annotated,subset=celltype %in% c("Neutrophil"))
DimPlot(selected_Neutrophil,reduction="umap",pt.size=1,label=T,split.by="group")
selected_Neutrophil.markers<-FindAllMarkers(selected_Neutrophil,only.pos=TRUE)
selected_Neutrophil.markers%>%
            group_by(cluster)%>%
            dplyr::filter(avg_log2FC>1)%>%
            slice_head(n=10)%>%
            ungroup()->top50       
Neutrophil_marker<-DoHeatmap(selected_Neutrophil,features=top50$gene)

#[18]:"Fam184b","Hsph1"
#[24]:"Crem","Maff","Nr4a3"
#[25]:"Itgae","Ptpn13","Pard3","Klf12"
Neutrophil_violin<-VlnPlot(selected_Neutrophil,features=c("Itgam","Sorl1","Lrp6","Itm2b","Igh-1a","H2-Ab1","Cd74","Abca1","Pkm",
"Nbr1","Notch"),split.by="group",ncol=3)#MHC-II gene Cd74 upregulated, Cited from <High-resolution analysis of the murine MHC class II> immunopeptidome
Neutrophil_violin<-VlnPlot(selected_Neutrophil,features=c("Ngp","Ltf","Pabpc1","Ifit3","Gm2a","Marco","Actg1","Mmp8","Apoa2","Ifit1","Spp1","Ccl4"),split.by="group",ncol=3)
ggsave("Neutrophil_violin.pdf",plot=Neutrophil_violin,width=24,height=20)

#######cDC#############
selected_cDC<-subset(Merged_Manually_Annotated,subset=celltype %in% c("cDC"))
DimPlot(selected_cDC,reduction="umap",pt.size=1,label=T,split.by="group",group.by="seurat_clusters")
selected_cDC.markers<-FindAllMarkers(selected_cDC,only.pos=TRUE)
selected_cDC.markers%>%
            group_by(cluster)%>%
            dplyr::filter(avg_log2FC>1)%>%
            slice_head(n=10)%>%
            ungroup()->top20
cDC_marker<-DoHeatmap(selected_cDC,features=top20$gene)
save(selected_cDC.markers,file="selected_cDC.markers.RData")
save(selected_cDC,file="selected_cDC.RData")
res<-gptcelltype(selected_cDC.markers,tissuename="human pbmc",model="gpt-3.5")
cDC_violin<-VlnPlot(selected_cDC,features=c("Pcdh7","Mfge8","Nbea","Sgcz","Nr3c2","Cadm2"),split.by="group",ncol=4)
#[4],[29]:"Pcdh7",Mfge8 ###adhesion and binding
#[9],[11],[14],[20],[22],[31]:"Npc2","Dstn","Clu"
#[13].[30]:"Gm42418","Lars2","Cdk8"
VlnPlot(selected_cDC,features=c("Gm42418","Lars2","Cdk8"))

NK_marker+pDC_marker+cDC_marker

ggsave("selected_markerheatmap.pdf",plot=NK_marker+ILC2_marker+cDC_marker,width=24,height=8)