## HER2-positive breast cancer subtypes
This repository contains functions and inputs required to compute HER2-positive subtypes from gene expression data as defined by Rediti M. et al. in PMID...
The R function to compute HER2-positive subtypes is available in "calc_HER2_groups_function.R".
The R script to compute HER2-positive subtypes in a simulated example gene expression matrix available in "example_compute_HER2_subtypes.R".

calc_HER2_groups() is the R function to compute the HER2-positive breast cancer subtypes, and requires the RDS objects "sigs_groups_class_final.RDS", "median_genes.RDS", "x_mean_genes.RDS", "x_sd_genes.RDS".

The scripts have been tested on the R software (v4.2.1).


### Information (details available in the .R scripts)
- Input is a gene expression matrix (gene symbols/entrezIDs as row names, sample IDs as column names), either RNA-seq TPM or FPKM normalized (not log-transformed) data, or microarray data already normalized (log2 transformed). 
It is possible to provide a list of gene expression datasets, which will be pre-processed separately with "median_rescale", then merged and further pre-processed with "standardize_data".
- The RDS files required (to be located in the home directory) are listed in the script and are available in this repository ("sigs_groups_class_final.RDS", "median_genes.RDS", "x_mean_genes.RDS", "x_sd_genes.RDS").
- The output is a dataframe with the score of each subtype, the subtype as category, and columns of each subtype vs rest to facilitate downstream analyses.

This classifier has not been tested with other normalization methods, as well as in datasets missing any of the genes used. Thus, while subtypes can still be identified, the accuracy in classifying samples may be reduced in these conditions.


### Dictionary of the subtypes:
- "IM" = Immune-enriched          
- "P/Met" = Proliferative/Metabolic-enriched  
- "Mes/S" = Mesenchymal/Stroma-enriched  
- "LUM" = Luminal  
- "ERBB2-E" = ERBB2-enriched  


### How to run the example:
The script to run the example is stored in "example_compute_HER2_subtypes.R" (details in the R script) and requires the content of "calc_HER2_groups_function.R", as well as "example_tpm_values.RDS" (an example gene expression matrix for testing, with hypothetical TPM-normalized gene expression data) and "HER2_subtypes_example_REF.RDS" (the expected results from the example).

Total runtime for the example gene expression matrix (genes included in the classifier for 500 samples) in "example_compute_HER2_subtypes.R" on a computer with 32GB RAM, Apple M1 Max CPU is ~0.03-0.06 seconds.

## Please cite our work when using these scripts and the present classification.

### Citation:




