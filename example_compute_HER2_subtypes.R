########### EXAMPLE SCRIPT

#### LOAD FILES required for CLASSIFIER: 
#### it requires "sigs_groups_class_final.RDS", "median_genes.RDS", "x_mean_genes.RDS", "x_sd_genes.RDS" in your working directory

sigs_groups_class_tot <- readRDS(paste0(getwd(), "/sigs_groups_class_final.RDS"))
median_genes = readRDS(paste0(getwd(), "/median_genes.RDS"))
xcenter_genes = readRDS(paste0(getwd(), "/x_mean_genes.RDS"))
xcenter_sd_genes = readRDS(paste0(getwd(), "/x_sd_genes.RDS"))

# Dictionary:
# "IM" = Immune-enriched         
# "P/Met" = Proliferative/Metabolic-enriched
# "Mes/S" = Mesenchymal/Stroma-enriched
# "LUM" = Luminal
# "ERBB2-E" = ERBB2-enriched
