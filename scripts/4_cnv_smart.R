#R
library(infercnv)
REF =c('GAMM','oligodendrocyte')
options(scipen=100)
obcnv <-CreateInfercnvObject(raw_counts_matrix='counts.matrix',
                             annotations_file='cells.txt', delim="\t",
                             ref_group_names=REF,
                             gene_order_file='feature.txt',
                             min_max_counts_per_cell=c(200, +Inf),
                             chr_exclude=c('chrX','chrY','chrM'))

#smart-seq2
options(error=function() traceback(2))
library(infercnv)
obcnv <-run(infercnv_obj=obcnv, cutoff=0.1,
            cluster_by_groups=TRUE,
            denoise=TRUE, num_threads=40,
            out_dir='CNVrun')
