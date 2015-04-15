# Null model for de novo regulatory mutations

# Probability of mutation in middle base of each trinucleotide is calculated from empirical data.
# The probabilities are assumed poisson, independent and summed across each sequence

mut_rates <- read.table("/Users/ps14/code/DDD/denovo_regs/data/forSanger_1KG_mutation_rate_table.txt", header=TRUE)

### SNP MODEL - based on sequence context ####

# the indel null model is inferred directly from global snp mutation rate and data - see rupit_core.R for generate_indel_null function

p_all <- function(from){
  
  # probability of mutation from base in position 2 to any of three other possible bases

  p = mut_rates$mu_snp[c(mut_rates$from == from)]
  return(sum(p))
}

p_sequence <- function(sequence){
  
  # sum p_all across each trinucleotide sliding along sequence
  
  sequence = as.character(sequence)
  p = lapply(seq(1:nchar(sequence)-2), function(z) p_all(substr(sequence, z, z+2)))
  return(sum(as.numeric(p)))
  
}


