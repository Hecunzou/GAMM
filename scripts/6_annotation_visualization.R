setwd('/public/home/hecun/glioma/GSE70630')
load('GSE70630_seurat.RData')
# add annotations
Idents(seura) <-'seurat_clusters'
levels(seura)
ides <-c("malignant_MGH54", "malignant_MGH53", "malignant_MGH36",
         "malignant_MGH54","malignant_MGH36","malignant_MGH60","GAMM",  
         "malignant_MGH60", "malignant_MGH60", "oligodendrocyte")

names(ides) <- levels(seura)
seura <-RenameIdents(seura, ides)

seura$celltype <-Idents(seura)


##markers for every cluster compared to all remaining cells, the positive cells
Idents(seura) <-'celltype'
allmarker <-FindAllMarkers(seura, only.pos=TRUE, test.use="wilcox",
                           min.pct=0.3, logfc.threshold=0.8)

write.csv(allmarker, file="GSE70630_allmarker.csv")



#visualize marker expression across clusters
feature="LINC01094"

#adjust contrast levels based on quantiles of non-zero expression
FeaturePlot(seura, features=feature,
            min.cutoff='q5', max.cutoff='q95',
            reduction="tsne", label=F)


#percentage of marker expression in each cluster
gene=c('LINC01094','LINC01736','NEAT1','FCGR1BP','FCGR1CP','HLA-DRB6','NAPSB','NCF1C',
      'C1QA','C3','CSF1R','GPR34','IFNGR1','IL1B','P2RY12','RGS10','TYROBP',
      'AIF1','CCL3','CCL4','CD14','CD53','CD68','CD163','CXCL8',
      'FCER1G','FCGR3A','HLA-DRA','VSIG4')

DotPlot(seura, features =unique(gene),
        group.by='celltype') + RotatedAxis()



