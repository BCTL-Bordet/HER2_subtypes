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


########### 
# First run functions stored in "calc_HER2_groups_function.R" 
###########

####### Installation of required packages

if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!require("matrixStats", quietly = TRUE))
  install.packages("matrixStats")


### LOAD THE EXAMPLE MATRIX INCLUDING THE GENES IN THE CLASSIFIER (SIMULATED TPM VALUES) and 500 SAMPLES, from you working directory
example_tpm_values <- readRDS(paste0(getwd(), "/example_tpm_values.RDS"))



## Example and duration
# start.time <- Sys.time()

library(dplyr)
library(matrixStats)

HER2_subtypes_example <- calc_HER2_groups(example_tpm_values, type = "TPM")

# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken

table(HER2_subtypes_example$HER2_subtype)


### compare output

## load the REFERENCE file to compare results, from your working directory
HER2_subtypes_example_REF <- readRDS(paste0(getwd(), "/HER2_subtypes_example_REF.RDS"))


head(HER2_subtypes_example_REF)
#             IM_score   P_Met_score Mes_S_score    LUM_score ERBB2_E_score HER2_subtype IM_vs_rest P_Met_vs_rest Mes_S_vs_rest LUM_vs_rest ERBB2_E_vs_rest
# Sample1 -0.002751829 -0.0081196634 -0.03497819  0.030171464   0.015678216          LUM          N             N             N           Y               N
# Sample2  0.030893550 -0.0688037377  0.03050687  0.018704724  -0.011301408           IM          Y             N             N           N               N
# Sample3 -0.003718349  0.0045502666  0.01431830 -0.011858745  -0.003291475        Mes/S          N             N             Y           N               N
# Sample4 -0.018050805  0.0307767329 -0.01073494 -0.005534733   0.003543749        P/Met          N             Y             N           N               N
# Sample5 -0.014120818 -0.0003796478  0.01044108 -0.002556931   0.006616317        Mes/S          N             N             Y           N               N
# Sample6  0.027334281 -0.0431264157  0.03907716 -0.004221809  -0.019063215        Mes/S          N             N             Y           N               N



table(HER2_subtypes_example_REF$HER2_subtype)
# IM   P/Met   Mes/S     LUM ERBB2-E 
# 83     162     137      25      93 

table(HER2_subtypes_example$HER2_subtype)

