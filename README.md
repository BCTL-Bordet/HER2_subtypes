### HER2-positive breast cancer subtypes
The function to compute HER2-positive subtypes is available in "calc_HER2_groups_function.R".
The script to compute HER2-positive subtypes in a simulated example gene expression matrix available in "example_compute_HER2_subtypes.R".

# calc_HER2_groups() compute the HER2-positive breast cancer subtypes as described in PMID...

The scripts have been tested on the R software (v4.2.1).

### Information (details available in the .R scripts)
- Input is a gene expression matrix (gene symbols/entrezIDs as row names, sample IDs as column names), either TPM or FPKM normalized (not log-transformed), or microarray data already normalized (log2 transformed). 
It is possible to provide a list of gene expression datasets, which will be pre-processed separately with "median_rescale", then merged and further pre-processed with "standardize_data".
- The RDS files required (to be located in the home directory) are listed in the script and are available in this repository ("sigs_groups_class_final.RDS", "median_genes.RDS", "x_mean_genes.RDS", "x_sd_genes.RDS").
- The output is a dataframe with the score of each subtype, the subtype as category, and columns of each subtype vs rest to facilitate downstream analyses.

This classifier has not been tested with other normalization methods, as well as in datasets missing any of the genes used. Thus, while subtypes can still be identified, the accuracy in classifying samples may be reduced in these conditions.

### Dictionary:
#"IM" = Immune-enriched         
#"P/Met" = Proliferative/Metabolic-enriched
#"Mes/S" = Mesenchymal/Stroma-enriched
#"LUM" = Luminal
#"ERBB2-E" = ERBB2-enriched


### Example:
The script to run the example is stored in /scripts/compute_HER2_subtypes.R and requires the content of /scripts/calc_HER2_groups_function.R, as well as "list_examples.RDS" and "tot_number_reads.RDS".

The files "measures.RDS" and "df_measures.RDS" represent the expected outputs created with script_for_example.R and can be checked to verify the results obtained with the given examples.

Total runtime for script_for_example.R on a computer with 32GB RAM, Apple M1 Max CPU is ~1.6-2 seconds.

Please cite our work when using this script.





The script to run the example is stored in /R/scripts and requires the content of /R/functions, "list_examples.RDS" and "tot_number_reads.RDS".

The files "measures.RDS" and "df_measures.RDS" represent the expected outputs created with script_for_example.R and can be checked to verify the results obtained with the given examples.

Total runtime for script_for_example.R on a computer with 32GB RAM, Apple M1 Max CPU is ~1.6-2 seconds.

Please cite our work when using this script.

C
