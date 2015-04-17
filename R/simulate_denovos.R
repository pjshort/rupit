# simulating assignment of de novo reads to regulatory regions

# rationale: file of regions will be used to filter the denovos. only the non-coding denovos
# which are contained in the genomic regions provided by REGIONS_PATH. These regions are presumed
# to be well-covered. Variants are filtered out which have a denovogear prior <0.1.

# requirements for REGIONS_PATH:
  # 1. $chr, $start, $stop
  # 2. $region_id (can be gene name for burden-style analysis, numerical id, or chr.start.stop)
# OPTIONAl:
  # 1. if $p_null is provided, no priors will need to be determined
  # 2. if $p_null not provided, $seq must be provided to calculate prior from mutation rates

source("/Users/ps14/code/DDD/denovo_regs/R/rupit/R/rupit_core.R")

# load list of filtered de novos - run rupit_core.R function "filter_denovos" if not yet filtered
DENOVOS_PATH = "/Users/ps14/code/DDD/denovo_regs/R/rupit/data/de_novo_filtered.txt"
REGIONS_PATH = "/Users/ps14/code/DDD/denovo_regs/R/rupit/data/regions_annotated.txt"
iterations = 1000

# load data
de_novo_filtered <- read.table(DENOVOS_PATH, sep="\t", header=TRUE)
well_covered_regions <- read.table(REGIONS_PATH, sep="\t", header=TRUE)

# infer parameters for estimation - will use p-null (probability of mutation) for regions if provided
if (any(!(c("region_id", "p_snp_null", "p_indel_null", "n_snp", "n_indel") %in% names(well_covered_regions)))){
  print("Need to run pre_process.R on regions, de novos first!")
}

snp_total <- sum(well_covered_regions$n_snp)

assign_de_novos <- function(regions, region_idx, snp_total){
  v = vector(length = nrow(regions))
  idx = sample(region_idx, size = snp_total, replace = TRUE, prob = regions$p_relative)
  counts = table(idx)
  idx = sort(unique(idx))
  v[idx] = counts
  return(v)
}

simulate_de_novos <- function(regions, snp_total, iterations){
  
  # populates a sparse matrix with counts of de novos in each region for each simulation
  
  regions$p_relative = regions$p_snp_null/sum(regions$p_snp_null)
  region_idx = seq(1:nrow(regions))
  
  # rows will be region and column will be count of de novos in each iteration
  sim_output = replicate(iterations, assign_de_novos(regions, region_idx, snp_total))
  return(sim_output)
}


plot_burden <- function(){
  
  # plots histograms of the number of genomic regions with 1 de novo, 2+ de novos
  
  results = simulate_de_novos(well_covered_regions, snp_total, iterations = iterations)
  
  single_hits = colSums(results == 1)
  two_or_more = colSums(results > 1)
  three_or_more = colSums(results > 2)
  
  # single de novos
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  h1 = hist(single_hits, xlab="# of regions with 1 denovo", main="Simulating De Novos in Well-Covered Regions", breaks=seq(min(single_hits) - 0.5, max(single_hits)+0.5, 1), col="cyan", xaxt="n")
  axis(side=1, at = h1$mids, labels=seq(min(single_hits),max(single_hits),1))
  observed_single = nrow(well_covered_regions[well_covered_regions$n_snp == 1,])
  abline(v=observed_single, col="black", lty=3, lw=5)
  
  # double de novos
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  h2 = hist(two_or_more, xlab="# of regions with 2 denovos", main="Simulating De Novos in Well-Covered Regions", breaks=seq(min(two_or_more)-0.5, max(two_or_more)+0.5, 1), col="cyan", xaxt="n")
  axis(side=1, at = h2$mids, labels=seq(min(two_or_more),max(two_or_more),1))
  observed_double = nrow(well_covered_regions[well_covered_regions$n_snp == 2,])
  abline(v=observed_double, col="black", lty=3, lw=5)
  
}

plot_burden()

