#### Load packages ####
library(readxl); library(dplyr)

#### Read data #### 
path <- "/Users/elisunorig/Desktop/All the Things/Grad School/Research/OMEI:SmarTeeth/OMEI/Rat Study/Data/"

feces_raw <- read_xlsx(paste0(path, "Mucosal paper_OTU reads transposed_0114.xlsx"), sheet = "Feces")
feces_factor <- as.factor(feces_raw$`Group name`)
feces_taxa <- as.matrix(dplyr::select(feces_raw,-c("Group name", "Group ID", "Sample name")))

plaque_raw <- read_xlsx(paste0(path, "Mucosal paper_OTU reads transposed_0114.xlsx"), sheet = "Plaque")
plaque_factor <- as.factor(plaque_raw$`Group name`)
plaque_taxa <- as.matrix(dplyr::select(plaque_raw,-c("Group name", "Group ID", "#NAME")))


saliva_raw <- read_xlsx(paste0(path, "Mucosal paper_OTU reads transposed_0114.xlsx"), sheet = "Saliva")
saliva_factor <- as.factor(saliva_raw$`Group name`)
saliva_taxa <- as.matrix(dplyr::select(saliva_raw,-c("Group name", "Group ID", "#NAME")))

#### Define helper functions ####
CompCts_overallTest <- function(factor_mat, taxa_mat, 
                                adj_method = c("holm", "hochberg", "hommel", "bonferroni", 
                                               "BH", "BY", "fdr", "none")){
  res_names <- c("Statistic", "pval", "signif", "pval_adj", "signif_adj")
  test1_res <- data.frame(matrix(nrow = ncol(taxa_mat), ncol = length(res_names)))
  colnames(test1_res) <- res_names; rownames(test1_res) <- colnames(taxa_mat)
  for(i in 1:ncol(taxa_mat)){
    res_i <- kruskal.test(taxa_mat[,i] ~ factor_mat) 
    test1_res$Statistic[i] <- res_i$statistic; test1_res$pval[i] <- res_i$p.value
  }
  test1_res$pval_adj <- stats::p.adjust(test1_res[,"pval"], method = adj_method)
  
  test1_res$signif <- gtools::stars.pval(test1_res$pval)
  test1_res$signif_adj <- gtools::stars.pval(test1_res$pval_adj)
  
  return(test1_res)
}

pairwise_OTU_tests <- function(taxa_mat, raw_mat, mbiome, quiet = T){
  res_mat_names <- c("Microbiome", "Taxa", "GroupA", "GroupB", "Stat.", "pval", "Signif.")
  res_mat <- data.frame(matrix(ncol = length(res_mat_names)))
  colnames(res_mat) <- res_mat_names
  
  grp_combos <- t(combn(unique(raw_mat$`Group name`), 2))
  
  for(i in 1:ncol(taxa_mat)){
    res_mat_i <- data.frame("Microbiome" = rep(mbiome, nrow(grp_combos)),
                            "Taxa" = colnames(taxa_mat)[i],
                            "GroupA" = grp_combos[,1],
                            "GroupB" = grp_combos[,2],
                            "Stat." = rep(NA, nrow(grp_combos)),
                            "pval" = rep(NA, nrow(grp_combos)),
                            "Signif." = rep(NA, nrow(grp_combos)))
    
    
    for(j in 1:nrow(grp_combos)){
      if(quiet){
        suppressWarnings(tmp_res <- wilcox.test(taxa_mat[which(raw_mat$`Group name`==res_mat_i$GroupA[j]),i],
                                                taxa_mat[which(raw_mat$`Group name`==res_mat_i$GroupB[j]),i]))
      } else {
        tmp_res <- wilcox.test(taxa_mat[which(raw_mat$`Group name`==res_mat_i$GroupA[j]),i],
                               taxa_mat[which(raw_mat$`Group name`==res_mat_i$GroupB[j]),i])
      }
      res_mat_i$Stat.[j] <- tmp_res$statistic
      res_mat_i$pval[j] <- tmp_res$p.value
    }
    res_mat_i$Signif. <- gtools::stars.pval(res_mat_i$pval)
    res_mat <- rbind(res_mat, res_mat_i)
  }
  
  res_mat$pval_adj <- stats::p.adjust(res_mat$pval, method = "bonf")
  res_mat$Signif.adj. <- gtools::stars.pval(res_mat$pval_adj)
  return(res_mat[-1,]) # drop the first empty row
}


#### Fecal overall test ####
feces_test <- CompCts_overallTest(factor_mat = feces_factor, taxa_mat = feces_taxa, adj_method = "bonferroni")
knitr::kable(feces_test, digits = 4, col.names = c("Stat.", "p", "Signif.", "p_adj", "Signif. (adj.)"), align = "c",
             caption = "Kruskal-Wallis test for significant difference between treatment groups (df = 4) in fecal microbiome")


#### Plaque overall test ####
plaque_test <- CompCts_overallTest(factor_mat = plaque_factor, taxa_mat = plaque_taxa, adj_method = "bonferroni")
knitr::kable(plaque_test, digits = 4, col.names = c("Stat.", "p", "Signif.", "p_adj", "Signif. (adj.)"), align = "c",
             caption = "Kruskal-Wallis test for significant difference between treatment groups (df = 4) in plaque microbiome")


#### Salivary overall test ####
saliva_test <- CompCts_overallTest(factor_mat = saliva_factor, taxa_mat = saliva_taxa, adj_method = "bonferroni")
knitr::kable(saliva_test, digits = 4, col.names = c("Stat.", "p", "Signif.", "p_adj", "Signif. (adj.)"), align = "c",
             caption = "Kruskal-Wallis test for significant difference between treatment groups (df = 4) in salivary microbiome")


#### Fecal pairwise test ####
feces_pairwise <- pairwise_OTU_tests(taxa_mat = feces_taxa, raw_mat = feces_raw, mbiome = "fecal")

for(i in 1:length(unique(feces_pairwise$Taxa))){
  which.taxa <- unique(feces_pairwise$Taxa)[i]
  print(knitr::kable(feces_pairwise[which(feces_pairwise$Taxa == which.taxa),-c(1,2)], digits = 4,
                     col.names = c("Group A", "Group B", "Stat. (W)", "p", "Signif.", "p_adj", "Signif. (adj.)"),
                     row.names = F, align = "c", 
                     caption = paste("Pairwise Wilcoxon rank sum tests for difference in OTU counts of ",
                                     which.taxa, "in fecal microbiome")))
}

#### Plaque pairwise test ####
plaque_pairwise <- pairwise_OTU_tests(taxa_mat = plaque_taxa, raw_mat = plaque_raw, mbiome = "fecal")

for(i in 1:length(unique(plaque_pairwise$Taxa))){
  which.taxa <- unique(plaque_pairwise$Taxa)[i]
  print(knitr::kable(plaque_pairwise[which(plaque_pairwise$Taxa == which.taxa),-c(1,2)], digits = 4,
                     col.names = c("Group A", "Group B", "Stat. (W)", "p", "Signif.", "p_adj", "Signif. (adj.)"),
                     row.names = F, align = "c", 
                     caption = paste("Pairwise Wilcoxon rank sum tests for difference in OTU counts of ",
                                     which.taxa, "in plaque microbiome")))
}

#### Salivary pairwise test ####

saliva_pairwise <- pairwise_OTU_tests(taxa_mat = saliva_taxa, raw_mat = saliva_raw, mbiome = "fecal")

for(i in 1:length(unique(saliva_pairwise$Taxa))){
  which.taxa <- unique(saliva_pairwise$Taxa)[i]
  print(knitr::kable(saliva_pairwise[which(saliva_pairwise$Taxa == which.taxa),-c(1,2)], digits = 4,
                     col.names = c("Group A", "Group B", "Stat. (W)", "p", "Signif.", "p_adj", "Signif. (adj.)"),
                     row.names = F, align = "c", 
                     caption = paste("Pairwise Wilcoxon rank sum tests for difference in OTU counts of ",
                                     which.taxa, "in salivary microbiome")))
}
