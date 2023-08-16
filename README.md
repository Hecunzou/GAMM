# the feature profile of glioma-associated macrophages/microglia (GAMM) from multiple scRNA-seq data

Code associated with "The GAMM cell-population related non-coding feature drive the progression of gliomas"

## Preparation

* Download scRNA-seq data from [GEO](https://www.ncbi.nlm.nih.gov/geo) database. Move it to the (data/)(https://github.com/Hecunzou/GAMM/tree/master/data) directory (e.g. ```zcat file.gz | mv - data/.```). The samples information was summarized in Supplementary Tables.

* Download the scripts in the(scripts/)(https://github.com/Hecunzou/GAMM/tree/master/scripts) directory. Create a project or environment, and install the following packages (e.g. [Seurat](https://satijalab.org/seurat), [infercnv](https://github.com/broadinstitute/infercnv/wiki)). 

## Analysis of scRNA-seq data

* The snRNA-seq data were processed with Seurat package (v4.3.0) in R (v4.2.0). First, each gene detected in more than 50 cells was considered as expressed, and each cell with at least 200 expressed features were retained for creating Seurat object. 
* Viable single cells were identified based on total counts abundance, effective number of features, and suitable mitochondria fraction. 
* The count data were normalized using the NormalizeData function with LogNormalize method for each dataset. The highly variable genes (2000) were identified via FindVariableFeatures function with vst method. Data were transformed via ScaleData with regressed variable "percent.mito". 
* Principal component analysis was then performed with the highly variable features. Then, cell clusters were determined by a shared nearest neighbor (SNN) modularity optimization. The t-distributed Stochastic Neighbor Embedding (t-SNE) was applied to project single cells onto a two-dimensional map to discover heterogeneity among all cells. 
* Cell-type specific gene features were identified via the FindMarkers/FindAllMarkers function from Seurat package, using a wilcox test, the “min.pct” argument was set to 0.3, the “logfc.threshold” argument was set to 1.0. The expression distribution and percentage of feature was visualized by FeaturePlot and DotPlot.

## CNV inferrence
* Chromosomal CNVs (copy number variations) was inferred from scRNA-seq data using the infercnv tool (v1.16.0) on Linux server. 
* The genes were sorted by their genomic chromosomal regions, and chromosomal copy number of genes in each cell was estimated by moving average of the relative expression of the neighboring genes within each chromosome. 
* The “denoise” argument was set to TRUE, the “cutoff=0.1” worked well for 10X Genomics, and “cutoff=1.0” worked well for Smart-seq2. 
* Multiple cell-clusters with accumulated CNVs were considered as malignant cells. Simultaneously, the non-malignant cells without obvious CNVs were annotated to specific cell-types with canonical marker genes, and identified as GAMM, mature oligodendrocytes, T-cells, neurons and endothelial cells, and so on.


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.
