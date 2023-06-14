## modelling

This folder contains the scripts used to tune, train and test the models in the Servius et al. "Predicting Class Switch Recombination in B-cells from Antibody Repertoire Data" paper. 

#### Terminology used in filenames: 

* **CP**: COVID patient, referring to the models fit with the Hospitalised COVID-19 Patient dataset.
* **RSV**: Respiratory syncytial virus, referring to models fit with the Respiratory Syncytial Virus Live Challenge dataset.
* **L1DO**: Leave one donor out.

#### Executing notes:

* Parallel: These scripts are set to run in parallel. The number of cores used has been set in each script using the `numCores` variable which is the number of cores requested minus 1. This variable can be adjusted as needs require and set to one if it is not possible to parallelise on your machine. 
* Cross-Validation: These scripts are written to be stand alone and complete. It would be possible to separate the hyperparameter search from the training and testing if desired. Attention needs to be paid in order to ensure the seed used is the same across separated scripts.

| folder name | filename | Description | 
|:-------------|:----------|:-------------|
| `all-donor-cv` |  | Contains the scripts for the all donor (AD) cross-validation strategy combined with the AD test as detailed in the Servius et al. paper. |
|  | `CP_LR.R` | Logistic regression model training and fitting on CP.  |  
|  | `CP_RF-NSx.R` | Random forest with variable nodesize tuning, training and fitting on CP. |  
|  | `CP_RF.R` | Random forest tuning, training and fitting on CP. |
|  | `CP_SVM.R` | Support vector machine tuning, training and fitting on CP. |  
|  | `RSV-CP_LASSO.R` | LASSO logistic regression tuning, training and fitting on RSV and CP (separately). |  
|  | `RSV_LR.R` | Logistic regression model training and fitting on RSV. |  
|  | `RSV_RF-NSx.R` | Random forest with variable nodesize tuning, training and fitting on RSV. |  
|  | `RSV_RF.R` | Random forest tuning, training and fitting on RSV. |  
|  | `RSV_SVM.R` | Support vector machine tuning, training and fitting on RSV. |
| `leave-one-donor-out-cv` |   | Contains the scripts for the leave one donor out (L1DO) cross-validation strategy combined with the L1DO test as detailed in the Servius et al. paper. |
|  | `CP_L1DO_LR.R` | Logistic regression L1DO training and fitting on CP. |    
|  | `CP_L1DO_RF.R` | Random forest L1DO tuning, training and fitting on CP. |   
|  | `RSV-CP_L1DO_LASSO.R` | LASSO logistic regression L1DO tuning, training and fitting on RSV and CP. | 
|  | `RSV_L1DO_RF-NSx.R` | Random forest with variable nodesize L1DO tuning, training and fitting on RSV. |  
|  | `RSV_L1DO_SVM.R` | Support vector machine L1dO tuning, training and fitting on RSV. |
|  | `CP_L1DO_RF-NSx.R` | Random forest with variable nodesize L1DO tuning, training and fitting on CP. |  
|  | `CP_L1DO_SVM.R` | Support vector machine L1DO tuning, training and fitting on CP. |  
|  | `RSV_L1DO_LR.R` | Logistic regression L1DO training and fitting on RSV. |
|  | `RSV_L1DO_RF.R` | Random forest L1DO tuning, training and fitting on RSV. |




