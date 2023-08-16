#The scRNA-seq data GSE70630 in this study were downloaded from GEO databases
setwd('/public/home/hecun/glioma/GSE70630')
matrx <-read.table("GSE70630_OG_processed_data_v2.txt.gz", sep=',',
                   header=T, row.names=1, check.names=F)
dim(matrx)
#create Seurat object,
library(Seurat)
seura <-CreateSeuratObject(counts =matrx, assay ="RNA",
                           min.cells=30, min.features=200,
                           names.field=1, names.delim='_',
                           project='GSE70630')

dim(seura)

##QC
#filter cells based on 'nFeature_RNA' and 'nCount_RNA'
seura <-subset(seura, subset=nFeature_RNA >200 & nCount_RNA >500)

seura <-subset(seura, subset=nFeature_RNA <1.0e+4 & nCount_RNA <7.5e+4)

#get ERCC genes, add to metadata of Seurat object
seura[["percent.mito"]] <-PercentageFeatureSet(seura, pattern ="^MT-")

#get ribosomal genes, add to metadata
seura[["percent.ribo"]] <-PercentageFeatureSet(seura, pattern='^RP[SL][[:digit:]]')

#filter cells based on 'percent.mito' and 'percent.ribo'
seura <-subset(seura, subset=percent.mito <15.8 & percent.ribo <37)


save(seura, file='GSE70630_seurat.RData')



