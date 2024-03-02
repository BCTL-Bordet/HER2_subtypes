# HER2_subtypes
The function to compute HER2-positive subtypes is available in "calc_HER2_groups_function.R"
the script to compute HER2-positive subtypes in an example gene expression matrix available in 

calc_HER2_groups() compute the HER2-positive breast cancer subtypes as described in PMID...

Input is a gene expression matrix (gene symbols/entrezIDs as row names, sample IDs as column names), either TPM or FPKM normalized (not log-transformed), or microarray data already normalized (log2 transformed). 

It is possible to provide a list of gene expression datasets, which will be pre-processed separately with "median_rescale", then merged and further pre-processed with "standardize_data".
The RDS files required (to be located in the home directory) are listed in the script and are available in this repository ("sigs_groups_class_final.RDS", "median_genes.RDS", "x_mean_genes.RDS", "x_sd_genes.RDS").

The output is a dataframe with the score of each subtype, the subtype as category, and columns of each subtype vs rest to facilitate downstream analyses.

This classifier has not been tested with other normalization methods, as well as in datasets missing any of the genes used. Thus, while subtypes can still be identified, the accuracy in classifying samples may be reduced in these conditions.

Quick Example:
# gene_matrix is a gene expression matrix with values as TPM
df <- calc_HER2_groups(gene_matrix, type = "TPM")




Functions to compute BCR/TCR diversity measures from MiXCR outputs. The scripts have been tested on the R software (v4.0.5 and v4.2.1). The two simulated samples were used to produce Supplementary Figure 1 in Rediti M, et al. Immunological and clinicopathological features predict HER2-positive breast cancer prognosis in the neoadjuvant NeoALTTO and CALGB 40601 randomized trials. Nat Commun. 2023. doi: 10.1038/s41467-023-42635-2. PMID: 37923752.

First, download the files list_examples.RDS and tot_number_reads.RDS.

list_examples.RDS is a list with 2 simulated samples containing, in this case, only IG information, and it shows how MiXCR output should be formatted to run the script (explanations for the columns are available in /R/scripts/script_for_example.R)
tot_number_reads.RDS is dataframe with the total number of reads mapping to genes for the two samples, required for the normalization of BCR/TCR number of reads
The functions to compute BCR/TCR measures are stored in /R/functions.

The script to run the example is stored in /R/scripts and requires the content of /R/functions, "list_examples.RDS" and "tot_number_reads.RDS".

The files "measures.RDS" and "df_measures.RDS" represent the expected outputs created with script_for_example.R and can be checked to verify the results obtained with the given examples.

Total runtime for script_for_example.R on a computer with 32GB RAM, Apple M1 Max CPU is ~1.6-2 seconds.

Please cite our work when using this script.

C
