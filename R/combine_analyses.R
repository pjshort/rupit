# library of functions to combine different analyses e.g. phenotype and recurrent mutational burden

library(survcomp) # needed for combined p-value testing
ddg2p <- read.table("../data/DDG2P_freeze_with_gencode19_genomic_coordinates_20141118_fixed.txt", sep="\t", header=TRUE)

merge_pheno_geno_pvals <- function(enrichment_test, phenotype_test, by){
  
  # merge enrichment_test and phenotype_test by whatever was used to split (usually region_id or closest_gene)
  # adds a factor with p-val for fisher's combined test (and sort by this p-value)
  
  merged = merge(enrichment_test, phenotype_test, by = by)
  pvals = cbind(merged$p_snp_test, merged$hpo_similarity_p_value)
  merged$p_combined = apply(pvals, MARGIN = 1, FUN = combine.test)
  merged = merged[order(merged$p_combined),]
  merged$in_ddg2p = merged$closest_gene %in% ddg2p$gencode_gene_name
  return(merged)
}