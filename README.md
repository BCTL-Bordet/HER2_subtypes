# HER2-positive breast cancer molecular subtypes
This repository contains functions and inputs required to compute HER2-positive breast cancer molecular subtypes from gene expression data as defined by Rediti M. _et al_. in "Identification of HER2-positive breast cancer molecular subtypes with potential clinical implications in the ALTTO clinical trial", PMID: 39613746 (https://www.nature.com/articles/s41467-024-54621-3). 
**Please note that our molecular classification was developed, applied and tested in HER2-positive breast cancer as defined by immunohistochemistry and fluorescence in situ hybridization criteria, including both hormone receptor-positive and negative tumors. Therefore, please select cohorts of HER2-positive breast cancer according to such criteria before computing the molecular subtypes with our classifier.**

The R function **calc_HER2_groups()** to compute HER2-positive breast cancer subtypes is available in "**R/calc_HER2_groups_function.R**".  
The R script to compute HER2-positive subtypes in a simulated example gene expression matrix available in "example_files/example_compute_HER2_subtypes.R".

**calc_HER2_groups()** requires the RDS objects "**sigs_groups_class_final.RDS**", "**median_genes.RDS**", "**x_mean_genes.RDS**", "**x_sd_genes.RDS**" contained in the directory "**input_files/**".

The scripts have been tested on the R software (v4.2.1).


# Information (details also available in the .R scripts)
- The main input is a **gene expression matrix** (gene symbols/entrezIDs as row names, sample IDs as column names, which does not need to be filtered for the genes included in the classifier), either RNA-seq TPM or FPKM normalized (not log-transformed) data, or microarray data already normalized (log2 transformed). 
It is possible to provide a list of gene expression datasets, which will be pre-processed separately with "median_rescale", then merged and further pre-processed with "standardize_data".
- The **RDS files** required (to be located in your working directory, or with _path/to/file/_ to be specified directly within calc_HER2_groups()) are listed in the script and are available in this repository in the directory "**input_files/**" ("sigs_groups_class_final.RDS", "median_genes.RDS", "x_mean_genes.RDS", "x_sd_genes.RDS"). They can also be saved as variables in your environment and called replacing the .RDS files in calc_HER2_groups().
- With the RDS files in your working directory, you can run the function on TPM gene expression data simply with:
  
  **subtypes_df <- calc_HER2_groups(your_TPM_matrix, type = "TPM")**
- The **output** is a data frame with the sample names as row names, scores for all subtype in each sample (signatures computed as weighted mean of the genes, using the coefficients derived from the LASSO classifier as weights, and after pre-processing as described in the publication and in the R scripts), the subtype as category assigned to each sample (based on the highest score), and columns of each subtype (as category) vs. rest to facilitate downstream analyses.

This classifier has not been tested with other normalization methods, as well as in datasets missing any of the genes used. Thus, while subtypes can still be identified, the accuracy in classifying samples may be reduced in these conditions.


## Dictionary of the subtypes:
- "IM" = Immune-enriched          
- "P/Met" = Proliferative/Metabolic-enriched  
- "Mes/S" = Mesenchymal/Stroma-enriched  
- "LUM" = Luminal  
- "ERBB2-D" = ERBB2-dependent  


# How to run the example:
The script and files to run the example are stored in the directory "example_files": the R script "**example_compute_HER2_subtypes.R**" (details to run it also reported in the R script) requires the content of "**R/calc_HER2_groups_function.R**", the RDS files in "**input_files/**", as well as "**example_tpm_values.RDS**" (an example gene expression matrix for testing, with hypothetical TPM-normalized gene expression data) and "**HER2_subtypes_example_REF.RDS**" (the expected results from the example).

Total runtime for the example gene expression matrix (87 genes included in the classifier in rows, 500 samples in columns) in "example_compute_HER2_subtypes.R" on a computer with 32GB RAM, Apple M1 Max CPU is ~0.03-0.06 seconds.

## Please cite our work when using the present classification and these scripts.

### Citation:
https://pubmed.ncbi.nlm.nih.gov/39613746/

### Further information
The molecular classification and the classifier were created by Mattia Rediti (MD, PhD) at the Breast Cancer Translational research Laboratory (BCTL), Institut Jules Bordet, ULB, in collaboration with David Venet (PhD) and under the supervision of Christos Sotiriou (MD, PhD).




