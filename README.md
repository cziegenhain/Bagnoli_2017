# Bagnoli_2017
This repository contains the code for analysis of single-cell RNA-seq data of Bagnoli et al., 2017.

## mcSCRB-seq: sensitive and powerful single-cell RNA sequencing
Single-cell RNA sequencing (scRNA-seq) has emerged as the central genome-wide method to characterize cellular identities and processes. While performance of methods is improving, an optimum in terms of sensitivity, efficiency and/or flexibility has not been reached yet. Among the plate-based counting methods that are well suited for quantifying gene expression for hundreds of cells across many samples, “Single-Cell RNA-Barcoding and Sequencing” (SCRB-seq) is one of the most powerful ones. Based on this protocol, we systematically evaluated reverse transcriptases, buffer modifications and PCR polymerases.
In particular, the addition of polyethylene glycol increased the sensitivity considerably. Based on these evaluations, we developed molecular crowding SCRB-seq (mcSCRB-seq), a fast, cost-efficient and sensitive protocol. By analyzing mouse embryonic stem cells and ERCC spike-ins we show that mcSCRB-seq is the most sensitive scRNA-seq method to date.

## Preprocessing
All scRNA-seq data was preprocessed with [zUMIs](https://github.com/sdparekh/zUMIs/) (Parekh et al., 2017).

The command was run as follows:

```Shell
bash zUMIs-master.sh -f JM8.read1.fastq.gz -r JM8.read2.fastq.gz -n JM8 -g /data/ngs/genomes/Mouse/mm10/STAR5idx_ERCC_noGTF/ -a /data/ngs/genomes/Mouse/mm10/Mus_musculus.GRCm38.75.clean.spike.gtf -c 1-14 -m 15-24 -l 50 -z 2 -u 3 -p 16 -R yes -d 10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000,3000000,4000000,5000000
```

This resulted in the [JM8.rds object found in this repository](Data/JM8.rds).


## Session Info
```R
> sessionInfo()
R version 3.4.0 (2017-04-21)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS Sierra 10.12.6

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib

locale:
[1] de_DE.UTF-8/de_DE.UTF-8/de_DE.UTF-8/C/de_DE.UTF-8/de_DE.UTF-8

attached base packages:
 [1] stats4    grid      parallel  splines   stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] hexbin_1.27.1       biomaRt_2.32.0      RColorBrewer_1.1-2  scales_0.5.0        scran_1.4.4        
 [6] BiocParallel_1.10.1 scater_1.4.0        Biobase_2.36.2      BiocGenerics_0.22.0 edgeR_3.18.1       
[11] limma_3.32.2        matrixStats_0.52.2  ineq_0.2-13         Hmisc_4.0-3         Formula_1.2-1      
[16] survival_2.41-3     lattice_0.20-35     bbmle_1.0.19        powsimRDev_0.0.905  doMC_1.3.4         
[21] iterators_1.0.8     foreach_1.4.3       gamlss.dist_5.0-2   bindrcpp_0.2        MASS_7.3-47        
[26] cowplot_0.8.0       dplyr_0.7.2         ggplot2_2.2.1      

loaded via a namespace (and not attached):
  [1] SparseM_1.77               rtracklayer_1.36.3         ggthemes_3.4.0             R.methodsS3_1.7.1         
  [5] lavaan_0.5-23.1097         coda_0.19-1                nonnest2_0.4-1             tidyr_0.7.0               
  [9] acepack_1.4.1              bit64_0.9-7                knitr_1.16                 irlba_2.2.1               
 [13] aroma.light_3.6.0          DelayedArray_0.2.7         R.utils_2.5.0              Rook_1.1-1                
 [17] data.table_1.10.5          rpart_4.1-11               hwriter_1.3.2              RCurl_1.95-4.8            
 [21] doParallel_1.0.10          snow_0.4-2                 GenomicFeatures_1.28.2     RSQLite_2.0               
 [25] VGAM_1.0-4                 combinat_0.0-8             bit_1.1-12                 httpuv_1.3.5              
 [29] ggsci_2.7                  SummarizedExperiment_1.6.3 DrImpute_1.0               assertthat_0.2.0          
 [33] viridis_0.4.0              tximport_1.4.0             RMTstat_0.3                IHW_1.4.0                 
 [37] caTools_1.17.1             igraph_1.0.1               DBI_0.7                    geneplotter_1.54.0        
 [41] htmlwidgets_0.8            EDASeq_2.10.0              RcppArmadillo_0.7.960.1.2  purrr_0.2.3               
 [45] backports_1.1.0            DDRTree_0.1.5              pbivnorm_0.6.0             permute_0.9-4             
 [49] scDD_1.0.0                 annotate_1.54.0            moments_0.14               RcppParallel_4.3.20       
 [53] blockmodeling_0.1.8        Cairo_1.5-9                quantreg_5.33              abind_1.4-5               
 [57] withr_1.0.2                RcppEigen_0.3.3.3.0        checkmate_1.8.2            GenomicAlignments_1.12.1  
 [61] fdrtool_1.2.15             mclust_5.3                 SCnorm_0.99.7              mnormt_1.5-5              
 [65] cluster_2.0.6              DEDS_1.50.0                NBPSeq_0.3.0               lazyeval_0.2.0            
 [69] crayon_1.3.2               genefilter_1.58.1          glmnet_2.0-10              pkgconfig_2.0.1           
 [73] slam_0.1-40                labeling_0.3               GenomeInfoDb_1.12.1        nlme_3.1-131              
 [77] vipor_0.4.5                devtools_1.13.2            nnet_7.3-12                bindr_0.1                 
 [81] rlang_0.1.2                miniUI_0.1.1               MatrixModels_0.4-1         sandwich_2.3-4            
 [85] extRemes_2.0-8             BPSC_0.99.1                cidr_0.1.5                 distillery_1.0-2          
 [89] Matrix_1.2-9               BASiCS_0.7.30              lpsymphony_1.4.1           zoo_1.8-0                 
 [93] base64enc_0.1-3            beeswarm_0.2.3             pheatmap_1.0.8             viridisLite_0.2.0         
 [97] rjson_0.2.15               bitops_1.0-6               shinydashboard_0.6.0       NOISeq_2.20.0             
[101] R.oo_1.21.0                Lmoments_1.2-3             spam_1.4-0                 KernSmooth_2.23-15        
[105] ggExtra_0.7                Biostrings_2.44.1          EBSeq_1.16.0               blob_1.1.0                
[109] rgl_0.98.1                 stringr_1.2.0              qvalue_2.8.0               msir_1.3.1                
[113] brew_1.0-6                 arm_1.9-3                  ShortRead_1.34.0           NbClust_3.0               
[117] S4Vectors_0.14.3           memoise_1.1.0              magrittr_1.5               plyr_1.8.4                
[121] gplots_3.0.1               gdata_2.18.0               zlibbioc_1.22.0            compiler_3.4.0            
[125] HSMMSingleCell_0.110.0     pcaMethods_1.68.0          lme4_1.1-13                DESeq2_1.16.1             
[129] fitdistrplus_1.0-9         Rsamtools_1.28.0           ade4_1.7-8                 DSS_2.16.0                
[133] XVector_0.16.0             htmlTable_1.9              mgcv_1.8-17                ROTS_1.4.0                
[137] MAST_1.2.1                 stringi_1.1.5              densityClust_0.2.1         locfit_1.5-9.1            
[141] latticeExtra_0.6-28        tools_3.4.0                monocle_2.4.0              foreign_0.8-67            
[145] outliers_0.14              bsseq_1.12.1               gridExtra_2.2.1            Rtsne_0.13                
[149] digest_0.6.12              FNN_1.1                    shiny_1.0.5                qlcMatrix_0.9.5           
[153] quadprog_1.5-5             Rcpp_0.12.12               car_2.1-4                  GenomicRanges_1.28.3      
[157] pscl_1.4.9                 AnnotationDbi_1.38.2       minpack.lm_1.2-1           colorspace_1.3-2          
[161] XML_3.98-1.7               fields_9.0                 IRanges_2.10.2             statmod_1.4.29            
[165] flexmix_2.3-14             xtable_1.8-2               nloptr_1.0.4               jsonlite_1.5              
[169] baySeq_2.10.0              dynamicTreeCut_1.63-1      modeltools_0.2-21          testthat_1.0.2            
[173] R6_2.2.2                   clusterCrit_1.2.7          htmltools_0.3.6            mime_0.5                  
[177] minqa_1.2.4                glue_1.1.1                 DT_0.2                     DESeq_1.28.0              
[181] RUVSeq_1.10.0              codetools_0.2-15           maps_3.1.1                 mvtnorm_1.0-6             
[185] tibble_1.3.4               pbkrtest_0.4-7             numDeriv_2016.8-1          ggbeeswarm_0.5.3          
[189] scde_1.99.4                gtools_3.5.0               openxlsx_4.0.17            CompQuadForm_1.4.3        
[193] cobs_1.3-3                 fastICA_1.2-0              munsell_0.4.3              rhdf5_2.20.0              
[197] GenomeInfoDbData_0.99.0    reshape2_1.4.2             gtable_0.2.0               NBGOF_0.2.2    
```
