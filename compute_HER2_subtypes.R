# color palette
col_HER2_subtypes <- c("IM" = "#db6d00",         
                       "P/Met" = "#920000",
                       "Mes/S" = "#009292",
                       "LUM" = "#006ddb",
                       "ERBB2-E" = "#ff6db6")

# it requires "sigs_groups_class_final.RDS", "median_genes.RDS", "x_mean_genes.RDS", "x_sd_genes.RDS" in your home directory
# output is a dataframe with the score of each subtype, the subtype as category ("HER2_subtype"), and columns of each subtype vs rest to facilitate downstream analyses
calc_HER2_groups <- function(d,  # d = gene expression data (rows are gene symbols/entrezIDs, columns are samples)
                             sig = readRDS("sigs_groups_class_final.RDS"), # signature file (a list of signatures with gene symbols, entrezID and coefficients)
                             type = c("TPM", "FPKM", "microarray"), # to be specified: after ratio-based correction, if TPM/FPKM, it transforms in log2(d+1), if microarray, assumed in log scale, so first un-logged 2^d and then log2(d+1). If standardize_data = TRUE and TPM/FPKM, centering and scaling are performed using center_genes and sd_genes. If Microarray, scale() is applied with default options
                             is_list_gene_exp = FALSE, # if TRUE, it allows to apply the same median ratio procedure to a list of gene expression matrices (TPM/FPKM) independently, which are then merged at the end before standardization
                             median_genes = readRDS("median_genes.RDS"), # the file containing median levels of TPM normalized gene expression from ALTTO TPM for genes included in the signatures
                             median_rescale = TRUE, # ratio-based correction: it will rescale gene expression data so that the median of selected genes will be the same as the median of TPM normalized genes in ALTTO (inspired by PMID: 28610557)
                             standardize_data_type = TRUE,  # if TRUE, it will apply scale() on gene expression data based on type, using scale(...,center=center_genes,scale=sd_genes) if TPM/FPKM, or scale(...,center=TRUE,scale=TRUE) if microarray. If FALSE, it will apply scale(...,center=TRUE,scale=TRUE) regardless
                             center_genes = readRDS("x_mean_genes.RDS"),
                             sd_genes = readRDS("x_sd_genes.RDS"),
                             fpkm_to_tpm = TRUE # it converts FPKM to TPM if type = "FPKM"
                             ) 
{ 
  '%!in%' <- function(x,y)!('%in%'(x,y))
  
  if(!require(dplyr)) { stop("dplyr library not installed"); }
  if (is.list(sig) && !is.data.frame(sig))
  { ret = sapply(sig, function(i) calc_HER2_groups(d, i, 
                                                   type = type, 
                                                   is_list_gene_exp = is_list_gene_exp, 
                                                   median_rescale = median_rescale, 
                                                   median_genes = median_genes,
                                                   standardize_data_type = standardize_data_type,  
                                                   center_genes = center_genes,
                                                   sd_genes = sd_genes, 
                                                   fpkm_to_tpm = fpkm_to_tpm));
  
  
  ret <- as.data.frame(ret)
  factor_tot_max <- as.data.frame(colnames(ret)[apply(ret,1,which.max)])
  
  factor_tot_max <- as.data.frame(factor_tot_max)
  colnames(factor_tot_max) <- "HER2_subtype"
  
  factor_tot_max$HER2_subtype <- factor(factor_tot_max$HER2_subtype, levels = c("IM", "P/Met", "Mes/S", "LUM", "ERBB2-E"))
  
  ret <- cbind(ret, factor_tot_max)
  colnames(ret) <- c("IM_score", "P_Met_score", "Mes_S_score", "LUM_score", "ERBB2_E_score", "HER2_subtype")
  ret$IM_vs_rest <- case_when(ret$HER2_subtype %in% "IM" ~ "Y",
                              ret$HER2_subtype %!in% "IM" ~ "N")
  ret$P_Met_vs_rest <- case_when(ret$HER2_subtype %in% "P/Met" ~ "Y",
                                 ret$HER2_subtype %!in% "P/Met" ~ "N")
  ret$Mes_S_vs_rest <- case_when(ret$HER2_subtype %in% "Mes/S" ~ "Y",
                                 ret$HER2_subtype %!in% "Mes/S" ~ "N")
  ret$LUM_vs_rest <- case_when(ret$HER2_subtype %in% "LUM" ~ "Y",
                               ret$HER2_subtype %!in% "LUM" ~ "N")
  ret$ERBB2_E_vs_rest <- case_when(ret$HER2_subtype %in% "ERBB2-E" ~ "Y",
                                   ret$HER2_subtype %!in% "ERBB2-E" ~ "N")
  
  return(ret);
  }
  
  if(median_rescale){
    if(type == "TPM"){
      if(is_list_gene_exp == FALSE){
        
        d_median <- colMedians(t(d), na.rm = TRUE)
        names(d_median) <- rownames(d)
        
        d_median <- d_median[!d_median == 0] # filter genes with median of 0
        
        d <- d[rownames(d) %in% names(d_median), ]
        
        d <- d[intersect(rownames(d), names(median_genes)), ]
        
        # ensure that genes are filtered if present and ordered
        median_genes_in_d <- median_genes[names(median_genes) %in% rownames(d)] # ALTTO genes in dataset
        d_median_in_altto <- d_median[names(d_median) %in% names(median_genes_in_d)] # Dataset genes in ALTTO
        median_genes_in_d <- median_genes_in_d[names(d_median_in_altto)]
        d_median_in_altto <- d_median_in_altto[names(median_genes_in_d)]
        
        d <- d[rownames(d) %in% names(d_median_in_altto), ] # dataset genes in common
        d <- d[names(d_median_in_altto), ]
        
        # compute ratio of medians
        median_genes_RATIO <- median_genes_in_d/d_median_in_altto
        
        df_median <- t(mapply("*", as.data.frame(t(d)), median_genes_RATIO))
        
        df_median <- as.matrix(df_median)
        
        colnames(df_median) <- colnames(d)
        
        d <- log2(df_median + 1)
        
      }
      
      if(is_list_gene_exp == TRUE){
        results <- list()
        for(i in 1:length(d)){
          
          d_list <- data.matrix(d[[i]])
          if (!is.matrix(d_list)) {
            stop("Error: Input data is not a matrix.")
          }
          
          d_median <- colMedians(t(d_list), na.rm = TRUE)
          names(d_median) <- rownames(d_list)
          
          d_list <- d_list[names(d_median), ]
          
          d_list <- d_list[intersect(rownames(d_list), names(median_genes)), ]
          
          # ensure that genes are filtered if present and ordered
          median_genes_in_d <- median_genes[names(median_genes) %in% rownames(d_list)] # ALTTO genes in dataset
          d_median_in_altto <- d_median[names(d_median) %in% names(median_genes_in_d)] # Dataset genes in ALTTO
          median_genes_in_d <- median_genes_in_d[names(d_median_in_altto)]
          d_median_in_altto <- d_median_in_altto[names(median_genes_in_d)]
          
          d_list <- d_list[rownames(d_list) %in% names(d_median_in_altto), ] # dataset genes in common
          d_list <- d_list[names(d_median_in_altto), ]
          
          # compute ratio of medians
          median_genes_RATIO <- median_genes_in_d/d_median_in_altto
          
          df_median <- t(mapply("*", as.data.frame(t(d_list)), median_genes_RATIO))
          
          df_median <- as.matrix(df_median)
          
          colnames(df_median) <- colnames(d_list)
          
          d_list <- log2(df_median + 1)
          results[[i]] <- d_list
        }
        
        # if matrices have different number of rows
        common_genes <- Reduce(intersect, lapply(results, rownames))
        # Genes not present in some of the gene expression matrices are filtered out
        filter_rows <- function(expression_matrix, common_genes) {
          expression_matrix <- expression_matrix[common_genes, , drop = FALSE]
          return(expression_matrix)
        }
        # Apply the function to each gene expression matrix
        results <- lapply(results, filter_rows, common_genes = common_genes)
        
        d <- do.call(cbind, results)
        
      } 
    }
    
    else if(type == "FPKM" ){
      
      if(fpkm_to_tpm){
      fpkmToTpm <- function(fpkm)
      {
        (fpkm / sum(fpkm)) * 10^6
      }
      
      fpkmToTpm_matrix <- function(fpkm_matrix) {
        fpkm_matrix_new <- apply(fpkm_matrix, 2, fpkmToTpm)  
      }
      
      if(is_list_gene_exp == FALSE){
        
        d <- fpkmToTpm_matrix(d)
        
        d_median <- colMedians(t(d), na.rm = TRUE)
        names(d_median) <- rownames(d)
        
        d_median <- d_median[!d_median == 0] # filter genes with median of 0
        
        d <- d[rownames(d) %in% names(d_median), ]
        
        d <- d[intersect(rownames(d), names(median_genes)), ]
        
        # ensure that genes are filtered if present and ordered
        median_genes_in_d <- median_genes[names(median_genes) %in% rownames(d)] # ALTTO genes in dataset
        d_median_in_altto <- d_median[names(d_median) %in% names(median_genes_in_d)] # Dataset genes in ALTTO
        median_genes_in_d <- median_genes_in_d[names(d_median_in_altto)]
        d_median_in_altto <- d_median_in_altto[names(median_genes_in_d)]
        
        d <- d[rownames(d) %in% names(d_median_in_altto), ] # dataset genes in common
        d <- d[names(d_median_in_altto), ]
        
        # compute ratio of medians
        median_genes_RATIO <- median_genes_in_d/d_median_in_altto
        
        df_median <- t(mapply("*", as.data.frame(t(d)), median_genes_RATIO))
        
        df_median <- as.matrix(df_median)
        
        colnames(df_median) <- colnames(d)
        
        d <- log2(df_median + 1)
        
      }
      
      if(is_list_gene_exp == TRUE){
        results <- list()
        for(i in 1:length(d)){
         
           d[[i]] <- fpkmToTpm_matrix(d[[i]])
          
          d_list <- data.matrix(d[[i]])
          if (!is.matrix(d_list)) {
            stop("Error: Input data is not a matrix.")
          }
          
          d_median <- colMedians(t(d_list), na.rm = TRUE)
          names(d_median) <- rownames(d_list)
          
          d_list <- d_list[names(d_median), ]
          
          d_list <- d_list[intersect(rownames(d_list), names(median_genes)), ]
          
          # ensure that genes are filtered if present and ordered
          median_genes_in_d <- median_genes[names(median_genes) %in% rownames(d_list)] # ALTTO genes in dataset
          d_median_in_altto <- d_median[names(d_median) %in% names(median_genes_in_d)] # Dataset genes in ALTTO
          median_genes_in_d <- median_genes_in_d[names(d_median_in_altto)]
          d_median_in_altto <- d_median_in_altto[names(median_genes_in_d)]
          
          d_list <- d_list[rownames(d_list) %in% names(d_median_in_altto), ] # dataset genes in common
          d_list <- d_list[names(d_median_in_altto), ]
          
          # compute ratio of medians
          median_genes_RATIO <- median_genes_in_d/d_median_in_altto
          
          df_median <- t(mapply("*", as.data.frame(t(d_list)), median_genes_RATIO))
          
          df_median <- as.matrix(df_median)
          
          colnames(df_median) <- colnames(d_list)
          
          d_list <- log2(df_median + 1)
          results[[i]] <- d_list
        }
        # if matrices have different number of rows
        common_genes <- Reduce(intersect, lapply(results, rownames))
        # Genes not present in some of the gene expression matrices are filtered out
        filter_rows <- function(expression_matrix, common_genes) {
          expression_matrix <- expression_matrix[common_genes, , drop = FALSE]
          return(expression_matrix)
        }
        # Apply the function to each gene expression matrix
        results <- lapply(results, filter_rows, common_genes = common_genes)
        
        d <- do.call(cbind, results)
        
      } 
      }
      else { # option to not transform FPKM to TPM
        if(is_list_gene_exp == FALSE){
          
          d_median <- colMedians(t(d), na.rm = TRUE)
          names(d_median) <- rownames(d)
          
          d_median <- d_median[!d_median == 0] # filter genes with median of 0
          
          d <- d[rownames(d) %in% names(d_median), ]
          
          d <- d[intersect(rownames(d), names(median_genes)), ]
          
          # ensure that genes are filtered if present and ordered
          median_genes_in_d <- median_genes[names(median_genes) %in% rownames(d)] # ALTTO genes in dataset
          d_median_in_altto <- d_median[names(d_median) %in% names(median_genes_in_d)] # Dataset genes in ALTTO
          median_genes_in_d <- median_genes_in_d[names(d_median_in_altto)]
          d_median_in_altto <- d_median_in_altto[names(median_genes_in_d)]
          
          d <- d[rownames(d) %in% names(d_median_in_altto), ] # dataset genes in common
          d <- d[names(d_median_in_altto), ]
          
          # compute ratio of medians
          median_genes_RATIO <- median_genes_in_d/d_median_in_altto
          
          df_median <- t(mapply("*", as.data.frame(t(d)), median_genes_RATIO))
          
          df_median <- as.matrix(df_median)
          
          colnames(df_median) <- colnames(d)
          
          d <- log2(df_median + 1)
          
        }
        
        if(is_list_gene_exp == TRUE){
          results <- list()
          for(i in 1:length(d)){
            
            d_list <- data.matrix(d[[i]])
            if (!is.matrix(d_list)) {
              stop("Error: Input data is not a matrix.")
            }
            
            d_median <- colMedians(t(d_list), na.rm = TRUE)
            names(d_median) <- rownames(d_list)
            
            d_list <- d_list[names(d_median), ]
            
            d_list <- d_list[intersect(rownames(d_list), names(median_genes)), ]
            
            # ensure that genes are filtered if present and ordered
            median_genes_in_d <- median_genes[names(median_genes) %in% rownames(d_list)] # ALTTO genes in dataset
            d_median_in_altto <- d_median[names(d_median) %in% names(median_genes_in_d)] # Dataset genes in ALTTO
            median_genes_in_d <- median_genes_in_d[names(d_median_in_altto)]
            d_median_in_altto <- d_median_in_altto[names(median_genes_in_d)]
            
            d_list <- d_list[rownames(d_list) %in% names(d_median_in_altto), ] # dataset genes in common
            d_list <- d_list[names(d_median_in_altto), ]
            
            # compute ratio of medians
            median_genes_RATIO <- median_genes_in_d/d_median_in_altto
            
            df_median <- t(mapply("*", as.data.frame(t(d_list)), median_genes_RATIO))
            
            df_median <- as.matrix(df_median)
            
            colnames(df_median) <- colnames(d_list)
            
            d_list <- log2(df_median + 1)
            results[[i]] <- d_list
          }
          # if matrices have different number of rows
          common_genes <- Reduce(intersect, lapply(results, rownames))
          # Genes not present in some of the gene expression matrices are filtered out
          filter_rows <- function(expression_matrix, common_genes) {
            expression_matrix <- expression_matrix[common_genes, , drop = FALSE]
            return(expression_matrix)
          }
          # Apply the function to each gene expression matrix
          results <- lapply(results, filter_rows, common_genes = common_genes)
          
          d <- do.call(cbind, results)
          
          
        } 
      }
    }
    else if(type == "microarray") {
      if(is_list_gene_exp == FALSE){
      d <- 2^d # inverse log transformation
      
      d_median <- colMedians(t(d), na.rm = TRUE)
      names(d_median) <- rownames(d)
      
      d_median <- d_median[!d_median == 0] # filter genes with median of 0
      
      d <- d[rownames(d) %in% names(d_median), ]
      
      d <- d[intersect(rownames(d), names(median_genes)), ]
      
      # ensure that genes are filtered if present and ordered
      median_genes_in_d <- median_genes[names(median_genes) %in% rownames(d)] # ALTTO genes in dataset
      d_median_in_altto <- d_median[names(d_median) %in% names(median_genes_in_d)] # Dataset genes in ALTTO
      median_genes_in_d <- median_genes_in_d[names(d_median_in_altto)]
      d_median_in_altto <- d_median_in_altto[names(median_genes_in_d)]
      
      d <- d[rownames(d) %in% names(d_median_in_altto), ] # dataset genes in common
      d <- d[names(d_median_in_altto), ]
      
      # compute ratio of medians
      median_genes_RATIO <- median_genes_in_d/d_median_in_altto
      
      df_median <- t(mapply("*", as.data.frame(t(d)), median_genes_RATIO))
      
      df_median <- as.matrix(df_median)
      
      colnames(df_median) <- colnames(d)
      
      d <- log2(df_median + 1)
      }
      if(is_list_gene_exp == TRUE){
        results <- list()
        for(i in 1:length(d)){
          
          d_list <- data.matrix(d[[i]])
          if (!is.matrix(d_list)) {
            stop("Error: Input data is not a matrix.")
          }
          d_list <- 2^d_list # inverse log transformation
          
          d_median <- colMedians(t(d_list), na.rm = TRUE)
          names(d_median) <- rownames(d_list)
          
          d_median <- d_median[!d_median == 0] # filter genes with median of 0
          
          d_list <- d_list[names(d_median), ]
          
          d_list <- d_list[intersect(rownames(d_list), names(median_genes)), ]
          
          # ensure that genes are filtered if present and ordered
          median_genes_in_d <- median_genes[names(median_genes) %in% rownames(d_list)] # ALTTO genes in dataset
          d_median_in_altto <- d_median[names(d_median) %in% names(median_genes_in_d)] # Dataset genes in ALTTO
          median_genes_in_d <- median_genes_in_d[names(d_median_in_altto)]
          d_median_in_altto <- d_median_in_altto[names(median_genes_in_d)]
          
          d_list <- d_list[rownames(d_list) %in% names(d_median_in_altto), ] # dataset genes in common
          d_list <- d_list[names(d_median_in_altto), ]
          
          # compute ratio of medians
          median_genes_RATIO <- median_genes_in_d/d_median_in_altto
          
          df_median <- t(mapply("*", as.data.frame(t(d_list)), median_genes_RATIO))
          
          df_median <- as.matrix(df_median)
          
          colnames(df_median) <- colnames(d_list)
          
          d_list <- log2(df_median + 1)
          
          results[[i]] <- d_list
        }
        # if matrices have different number of rows
        common_genes <- Reduce(intersect, lapply(results, rownames))
        # Genes not present in some of the gene expression matrices are filtered out
        filter_rows <- function(expression_matrix, common_genes) {
          expression_matrix <- expression_matrix[common_genes, , drop = FALSE]
          return(expression_matrix)
        }
        # Apply the function to each gene expression matrix
        results <- lapply(results, filter_rows, common_genes = common_genes)
        
        d <- do.call(cbind, results)
    }
  }
  }
  
  else {
    d <- log2(d + 1)
  }
  
  if(standardize_data_type) {
    if(type == "FPKM" ){
      d <- d[names(center_genes)[names(center_genes) %in% rownames(d)], ]
      d <- d[match(rownames(d), names(center_genes)[names(center_genes) %in% rownames(d)]), ]
      d = t(scale(t(d), center = center_genes[names(center_genes) %in% rownames(d)], 
                  scale = sd_genes[names(sd_genes) %in% rownames(d)]))
    }
    else if(type == "TPM"){
      d <- d[names(center_genes)[names(center_genes) %in% rownames(d)], ]
      d <- d[match(rownames(d), names(center_genes)[names(center_genes) %in% rownames(d)]), ]
      d = t(scale(t(d), center = center_genes[names(center_genes) %in% rownames(d)], 
                  scale = sd_genes[names(sd_genes) %in% rownames(d)]))
    }
    else if(type == "microarray"){
      d = t(scale(t(d), center = TRUE, scale = TRUE))
    }
  }
  else {
    d = t(scale(t(d), center = TRUE, scale = TRUE))
  }
  
  s = as.character(sig[,1]);
  s[s==""] = NA;
  x = intersect(rownames(d), s);
  x = x[!is.na(x)];
  if (length(x) == 0)
  { if (colnames(sig)[1] == "name") { other = "entrez"; } else { other = "name"; }
    s = as.character(sig[,other]);
    s[s==""] = NA;
    x = intersect(rownames(d), s);
    x = x[!is.na(x)];
    if (length(x) == 0)
    { warning("No gene in common");
      return(rep(NA, ncol(d)));
    }
  }
  
  sig = sig[match(x, s),];
  
  fCalc = function(d, x, coef) { colMeans(d[x,,drop=FALSE]*coef); }

  val = fCalc(d, x, sig[,"coefficient"]);
  return(val)
  
}
