# HER2_subtypes
script to compute HER2-positive subtypes available in "compute_HER2_subtypes.R"

calc_HER2_groups() compute the HER2-positive breast cancer subtypes as described in PMID...

Input is a gene expression matrix (genes as row names, sample IDs as column names), either TPM or FPKM normalized, or microarray data already normalized. It has not been tested with other normalizations, as well as in datasets missing any of the genes used.
It is possible to provide a list of gene expression datasets, which will be pre-processed separately with "median_rescale", then merged and further pre-processed with "standardize_data".
The RDS files required (to be located in the home directory) are listed in the script and are available in this repository.

The output is a dataframe with the score of each subtype, the subtype as category, and columns of each subtype vs rest to facilitate downstream analyses.

Example:

df <- calc_HER2_groups(gene_matrix, type = "TPM")


