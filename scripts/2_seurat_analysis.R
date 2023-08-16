setwd('/public/home/hecun/glioma/GSE70630')
load('GSE70630_seurat.RData')
#seurat analysis
library(Seurat)
seura <-NormalizeData(object=seura, scale.factor=10000,
                      normalization.method="LogNormalize")

seura <-FindVariableFeatures(seura, selection.method="vst",
                             nfeatures =2000)
# plot variable features with labels
top <-head(VariableFeatures(seura), 15)
plot <-VariableFeaturePlot(seura)

#linear transformation
seura <-ScaleData(object=seura, vars.to.regress="percent.mito")

#PCA
seura <-RunPCA(object=seura, assay="RNA", verbose=FALSE,
               features=VariableFeatures(seura))

print(seura[['pca']], dims =1:15, nfeatures=20)

DimPlot(object=seura, reduction='pca', pt.size=0.1, label=F)

# determine PCA scores
seura <-JackStraw(object=seura, num.replicate=100)
seura <-ScoreJackStraw(seura, dims=1:20)

JackStrawPlot(seura, dims=1:20)

ElbowPlot(object=seura)

#cluster cells, Louvain algorithm
seura <-FindNeighbors(seura, reduction="pca", dims=1:20)

seura <-FindClusters(seura, algorithm=1, resolution=0.5)


#UMAP
seura <-RunUMAP(object=seura, dims=1:20)

DimPlot(seura, reduction ="umap", pt.size=0.1)

#tSNE
seura <-RunTSNE(object=seura, dims=1:20)

DimPlot(seura, reduction='tsne', pt.size=0.1)

save(seura, file='GSE70630_seurat.RData')



