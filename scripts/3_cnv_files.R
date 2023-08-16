setwd('/public/home/hecun/glioma/GSE70630')
load("GSE70630_seurat.RData")
options(digits=3)
library(Seurat)
matrx <-GetAssayData(seura, slot='counts')
str(matrx)
matrx <-as.matrix(matrx)
write.table(matrx, file='counts.matrix', sep='\t', 
            quote=F, col.names=T, row.names=T)

#细胞类型文件
meta <-data.frame(barcode =rownames(seura@meta.data),
                  cluster =seura@meta.data$seurat_clusters,
                  celltype =seura@meta.data$celltype)
meta <-meta[row.names(seura@meta.data), c(1,3)]

write.table(meta, file='cells.txt',sep='\t',
            quote=F, col.names=F, row.names=F)

#基因位置排序文件
feature <-data.frame(name=rownames(seura),
                     chr=NA, start=NA, end=NA)
#merge gene
feature$chr <-DFano$chr[match(feature$name, DFano$gene)]
feature$start <-DFano$start[match(feature$name, DFano$gene)]
feature$end <-DFano$end[match(feature$name, DFano$gene)]

eve =c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
       'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18',
       'chr19','chr20','chr21','chr22')
feature[ ,2] <-factor(as.character(feature[ ,2]), levels=eve)

write.table(feature, file='feature.txt',sep='\t',
            quote=F, col.names=F, row.names=F)




