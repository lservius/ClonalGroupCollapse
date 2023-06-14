## data-processing

This folder contains the scripts used to process the data in Servius et al. "Predicting Class Switch Recombination in B-cells from Antibody Repertoire Data" paper. 

| filename | Description | 
|----------|-------------|
| `dataProcessing_CGC.R` | The script used to process the Hospitalised COVID-19 Patient and Respiratory Syncytial Virus Live Challenge ([Stewart et al. Frontiers 2022](https://doi.org/10.3389/fimmu.2022.807104)) data, including the data cleaning and the clonal group collapse (CGC). |
| `functions_dataProcessing.R` | Contained are the additional functions that are called in the `dataProcessing_CGC.R` script to process the data. Also includes the function used for oversampling the data for the training sets. The oversampling function is, however, included in the individual scripts for the actual model tuning and training so no requirement to source this script.  |
