Supplementary information / reproducible research files for the manuscript 
Title: "Predicting class switch recombination in B-cells from antibody repertoire data"

Authors: Servius, L.; Pigoli, P.; Ng, J.; Fraternali, F.

Code written by Servius, L.

The code was written and run in R with the following software versions:
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.6 LTS

Matrix products: default
BLAS:   /software/spackages_prod/apps/linux-ubuntu20.04-zen2/gcc-9.4.0/r-4.1.1-w2swqhtlwhxzg45ipiuh2fnj4zefy7yh/rlib/R/lib/libRblas.so
LAPACK: /software/spackages_prod/apps/linux-ubuntu20.04-zen2/gcc-9.4.0/r-4.1.1-w2swqhtlwhxzg45ipiuh2fnj4zefy7yh/rlib/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8   
 [7] LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] viridis_0.6.5        viridisLite_0.4.0    cowplot_1.1.1       
 [4] plyr_1.8.7           forcats_0.5.1        stringr_1.4.0       
 [7] dplyr_1.0.9          purrr_0.3.4          readr_2.1.2         
[10] tidyr_1.2.0          tibble_3.1.7         ggplot2_3.3.6       
[13] tidyverse_1.3.1      e1071_1.7-9          randomForest_4.7-1.1
[16] glmnet_4.1-4         Matrix_1.5-3         doParallel_1.0.17   
[19] iterators_1.0.14     foreach_1.5.2       

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.8.3     lubridate_1.8.0  lattice_0.20-45  class_7.3-20    
 [5] assertthat_0.2.1 utf8_1.2.2       R6_2.5.1         cellranger_1.1.0
 [9] backports_1.4.1  reprex_2.0.1     httr_1.4.3       pillar_1.7.0    
[13] rlang_1.0.2      readxl_1.4.0     rstudioapi_0.13  splines_4.1.1   
[17] munsell_0.5.0    proxy_0.4-26     broom_0.8.0      compiler_4.1.1  
[21] modelr_0.1.8     pkgconfig_2.0.3  shape_1.4.6      tidyselect_1.1.2
[25] gridExtra_2.3    codetools_0.2-18 fansi_1.0.3      crayon_1.5.1    
[29] tzdb_0.3.0       dbplyr_2.1.1     withr_2.5.0      grid_4.1.1      
[33] jsonlite_1.8.0   gtable_0.3.0     lifecycle_1.0.1  DBI_1.1.2       
[37] magrittr_2.0.3   scales_1.2.0     cli_3.3.0        stringi_1.7.6   
[41] fs_1.5.2         xml2_1.3.3       ellipsis_0.3.2   generics_0.1.2  
[45] vctrs_0.4.1      tools_4.1.1      glue_1.6.2       hms_1.1.1       
[49] survival_3.3-1   colorspace_2.0-3 rvest_1.0.2      haven_2.5.0  

----

Please note the following acronyms used:

LR: logistic regression
RF: random forest
RF-NSx: random forest with variable nodesize
SVM: support vector machine
GLMM: Generalised Linear Mixed Effects Model

CP: Hospitalised COVID-19 Patient
RSV: Respiratory Syncytial Virus Live Challenge

----

This folder contains the scripts for the Servius et al. "Predicting Class Switch Recombination in B-cells from Antibody Repertoire Data" manuscript. 

.	main.py		This file produces the figures and tables from the paper by pulling from the ./results folder. Please note that it runs some scripts within the ./scripts folder, check to see that the numCores variables in these scripts work for your system before running.

./data
	dataProcessing_CGC.R	script used to process the Hospitalised COVID-19 Patient and Respiratory Syncytial Virus Live Challenge (Stewart et al. Frontiers 2022) data, including the data cleaning and the clonal group collapse (CGC) - the dataset will need to be downloaded in to the top level of this data folder to use this script.
	func_dataProcessing.R	contained are the additional functions that are called in the dataProcessing_CGC.R script to process the data. Also includes the function used for oversampling the data for the training, however this has been separately placed in the ./functions folder and will be referenced from there for the model scripts.
	
	/processed
			covid-cg.txt		processed Hospitalised COVID-19 Patient data from Stewart et al. Frontiers (2022) using dataProcessing_CGC.R script.
			rsv-cg.txt		processed Respiratory Syncytial Virus Live Challenge data from Stewart et al. Frontiers (2022) using dataProcessing_CGC.R script.
			
/clean_noCGC		contains the intermediate step for the data-processing at the cleaning stage, before clonal group collapse - used to produce part of Table S1 in main.py.
	

./functions
	func_data-balance.R	function used for oversampling for the minority class which indicates a cross-class clonal group.
	func_LASSO-CV-AD.R	function used to perform the all-donor (AD) cross-validation of the hyperparameter (lambda) search for LASSO logistic regression (LR).
	func_LASSO-CV-L1DO	function used to perform the leave-one-donor-out (L1DO) cross-validation of the hyperparameter (lambda) search for LASSO logistic regression (LR).
	func_predictThreshold.R	function used to convert the probability prediction from LR and LASSO LR into binary class prediction.


./reduced-scripts		
	COVID-RSV_LASSO_RF_RF-NSx_SVM.Rmd	markdown file with the "reduced" scripts for RF, RF-NSx and SVM scripts, these have the hyperparameter search removed from the scripts and instead loads them from the ./hyperparam folder.
	/hyperparam				contains the hyperparameters selected for the RF, RF-NSx and SVM scripts to reduce computation time.

./scripts
	CP_LR.R 		LR model training and fitting on CP.
	CP_RF-NSx.R 		RF-NSx tuning, training and fitting on CP.
	CP_RF.R 		RF tuning, training and fitting on CP.
	CP_SVM.R 		SVM tuning, training and fitting on CP.
	RSV-CP_LASSO.R 		LASSO LR tuning, training and fitting on RSV and CP (separately).
	RSV_LR.R 		LR model training and fitting on RSV.
	RSV_RF-NSx.R 		RF-NSx tuning, training and fitting on RSV.
	RSV_RF.R 		RF tuning, training and fitting on RSV.
	RSV_SVM.R 		SVM tuning, training and fitting on RSV.
	CP_GLMM.R		GLMM tuning, training and fitting on CP.
	RSV_GLMM.R		GLMM tuning, training and fitting on RSV.

	CP_L1DO_LR.R 		LR L1DO training and fitting on CP.
	CP_L1DO_RF.R 		RF L1DO tuning, training and fitting on CP.
	RSV-CP_L1DO_LASSO.R 	LASSO LR L1DO tuning, training and fitting on RSV and CP.
	RSV_L1DO_RF-NSx.R 	RF-NSx L1DO tuning, training and fitting on RSV.
	RSV_L1DO_SVM.R 		SVM L1dO tuning, training and fitting on RSV.
	CP_L1DO_RF-NSx.R 	RF-NSx L1DO tuning, training and fitting on CP.
	CP_L1DO_SVM.R 		SVM L1DO tuning, training and fitting on CP.
	RSV_L1DO_LR.R 		LR L1DO training and fitting on RSV.
	RSV_L1DO_RF.R 		RF L1DO tuning, training and fitting on RSV.
	CP_L1DO_GLMM.R          GLMM L1DO tuning, training and fitting on CP.
        RSV_L1DO_GLMM.R         GLMM L1DO tuning, training and fitting on RSV.
	RSV_L1DO_RFx_modInt.R	RF-NSx L1DO single run on RSV for the model interpretation section in the paper.
	CP_L1DO_RFx_modInt.R    RF-NSx L1DO single run on CP for the model interpretation section in the paper.

./results	the outputs for the scripts that are needed to produce the figures and tables in the manuscipt will be/are saved here. 
	/figures_tables		the figures and tables produced from the data for the manuscript and support information are saved here.

----

Notes on running the scripts:
./scripts
	* Parallel: These scripts are set to run in parallel. The number of cores used has been set in each script using the numCores variable which is the number of cores requested minus 1. This variable can be adjusted as needs require and set to one if it is not possible to parallelise on your machine.
	* Cross-Validation: It would be possible to separate the hyperparameter search from the training and testing if desired. Attention needs to be paid in order to ensure the seed used is the same across separated scripts. 
	* The scripts source the functions in the ./functions folder as well as the data in ./data folder.
	* The scripts all assume the working directory is where they are saved. All scripts in the ./scripts folder will run from ./scripts, please keep this in mind and add a setwd() at the beginning of the scripts if necessary.
	* Please check the outputs for both LR and LASSO scripts, the CSV file may need to be reformatted for it to be readable but the main.py script that knits together the results.
	* Table S6 in the manuscript is the user time output by wrapping the scripts in system.time(), the base R package.
