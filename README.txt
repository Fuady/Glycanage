Title:
Association analysis of menopausal status and biological age measured using glycan age index from different platforms

Authors:
Angga M. Fuady* , Department of Biomedical Data Sciences, Leiden University Medical Center, The Netherlands
Said el Bouhaddani, Department of Data Science and Biostatistics, UMC Utrecht, div. Julius Centre, Huispost Str. 6.131, 3508 GA, Utrecht, The Netherlands
Hae-Won Uh, Department of Data Science and Biostatistics, UMC Utrecht, div. Julius Centre, Huispost Str. 6.131, 3508 GA, Utrecht, The Netherlands
Jeanine J. Houwing-Duistermaat, Department of Statistics, University of Leeds, United Kingdom

* The author that mainly responsible for writing the code
Email: A.M.Fuady@lumc.nl

> sessionInfo()
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows >= 8 x64 (build 9200)

Matrix products: default

locale:
[1] LC_COLLATE=English_Netherlands.1252  LC_CTYPE=English_Netherlands.1252   
[3] LC_MONETARY=English_Netherlands.1252 LC_NUMERIC=C                        
[5] LC_TIME=English_Netherlands.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] reshape_0.8.7         PO2PLS_0.9.0          mice_3.6.0            lattice_0.20-38      
 [5] mvtnorm_1.0-8         xtable_1.8-4          systemfit_1.1-22      lmtest_0.9-36        
 [9] zoo_1.8-3             car_3.0-4             carData_3.0-2         Matrix_1.2-17        
[13] msm_1.6.7             OmicsPLS_1.2.3        preprocessCore_1.42.0 stringr_1.4.0        
[17] gplots_3.0.1          magrittr_1.5          ggplot2_3.2.1         dplyr_0.8.3          
[21] readr_1.1.1          

loaded via a namespace (and not attached):
 [1] nlme_3.1-137              matrixStats_0.55.0        bitops_1.0-6             
 [4] bit64_0.9-7               tools_3.5.1               backports_1.1.5          
 [7] R6_2.4.0                  rpart_4.1-13              KernSmooth_2.23-16       
[10] mgcv_1.8-30               DBI_1.0.0                 lazyeval_0.2.2           
[13] BiocGenerics_0.26.0       colorspace_1.4-1          jomo_2.6-10              
[16] nnet_7.3-12               withr_2.1.2               tidyselect_0.2.5         
[19] bit_1.1-14                curl_4.2                  compiler_3.5.1           
[22] Biobase_2.40.0            expm_0.999-4              sandwich_2.5-1           
[25] labeling_0.3              caTools_1.17.1.2          scales_1.0.0             
[28] genefilter_1.64.0         digest_0.6.27             foreign_0.8-72           
[31] minqa_1.2.4               rio_0.5.16                pkgconfig_2.0.3          
[34] lme4_1.1-21               limma_3.36.2              rlang_0.4.1              
[37] readxl_1.3.1              RSQLite_2.1.2             rstudioapi_0.10          
[40] generics_0.0.2            BiocParallel_1.14.2       gtools_3.8.1             
[43] zip_2.0.4                 RCurl_1.95-4.12           Rcpp_1.0.2               
[46] munsell_0.5.0             S4Vectors_0.18.3          abind_1.4-5              
[49] lifecycle_0.1.0           stringi_1.4.3             yaml_2.2.0               
[52] MASS_7.3-50               plyr_1.8.4                blob_1.2.0               
[55] grid_3.5.1                parallel_3.5.1            gdata_2.18.0             
[58] forcats_0.4.0             mitml_0.3-7               crayon_1.3.4             
[61] haven_2.1.1               splines_3.5.1             annotate_1.60.0          
[64] hms_0.5.2                 zeallot_0.1.0             pillar_1.4.2             
[67] boot_1.3-20               reshape2_1.4.3            stats4_3.5.1             
[70] pan_1.6                   XML_3.98-1.20             glue_1.4.2               
[73] RcppArmadillo_0.9.800.1.0 data.table_1.12.6         vctrs_0.2.0              
[76] nloptr_1.2.1              cellranger_1.1.0          gtable_0.3.0             
[79] purrr_0.3.4               tidyr_1.0.0               assertthat_0.2.1         
[82] openxlsx_4.1.2            broom_0.5.2               survival_2.42-3          
[85] tibble_2.1.3              IRanges_2.14.10           memoise_1.1.0            
[88] AnnotationDbi_1.42.1      sva_3.30.0   
