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

source("/Users/ps14/code/DDD/denovo_regs/R/rupit_core.R")

# load list of filtered de novos - run rupit_core.R function "filter_denovos" if not yet filtered
DENOVOS_PATH = "/Users/ps14/code/DDD/denovo_regs/R/rupit/data/de_novo_filtered.txt"
REGIONS_PATH = "/Users/ps14/code/DDD/denovo_regs/R/rupit/data/regions_annotated.txt"
iterations = 1000

# load data
de_novo_filtered <- read.table(DENOVOS_PATH, sep="\t", header=TRUE)
well_covered_regions <- read.table(REGIONS_PATH, sep="\t", header=TRUE)

# infer parameters for estimation - will use p-null (probability of mutation) for regions if provided
if (!("p_null" %in% names(well_covered_regions))){
  well_covered_regions$p_null = generate_prior(well_covered_regions)
}


# count n_snps per region from fed data
if (!("n_snp" %in% names(well_covered_regions))){
  counts = mapply(count_denovos, well_covered_regions$chr, well_covered_regions$start, well_covered_regions$stop, MoreArgs = list(de_novo_filtered))
  
  well_covered_regions$n_snp = counts[1,]
  well_covered_regions$n_indel = counts[2,]
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
  
  regions$p_relative = regions$p_null/sum(regions$p_null)
  region_idx = seq(1:nrow(regions))
  
  # rows will be region and column will be count of de novos in each iteration
  sim_output = replicate(iterations, assign_de_novos(regions, region_idx, snp_total))  
  return(sim_output)
}


main <- function(){
  results = simulate_de_novos(well_covered_regions, snp_total, iterations = iterations)
  return(results)
}

plot_burden <- function(){
  
  # plots histograms of the number of genomic regions with 1 de novo, 
  # number with 2 de novos, and number with 3+ de novos
  
  results = main()
  
  single_hits = colSums(results == 1)
  double_hits = colSums(results == 2)
  three_or_more = colSums(results > 2)
  
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  h1 = hist(single_hits, xlab="# of regions with 1 denovo", main="Simulating De Novos in Well-Covered Regions", breaks=seq(min(single_hits) - 0.5, max(single_hits)+0.5, 1), col="cyan", xaxt="n")
  axis(side=1, at = h1$mids, labels=seq(min(single_hits),max(single_hits),1))
  observed_single = nrow(well_covered_regions[well_covered_regions$n_snp == 1,])
  abline(v=observed_single, col="black", lty=3, lw=5)
  
  # double de novos
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  h2 = hist(double_hits, xlab="# of regions with 2 denovos", main="Simulating De Novos in Well-Covered Regions", breaks=seq(min(double_hits)-0.5, max(double_hits)+0.5, 1), col="cyan", xaxt="n")
  axis(side=1, at = h2$mids, labels=seq(min(double_hits),max(double_hits),1))
  observed_double = nrow(well_covered_regions[well_covered_regions$n_snp == 2,])
  abline(v=observed_double, col="black", lty=3, lw=5)
  
  # triple de novos
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  h3 = hist(three_or_more, xlab="# of regions with 3+ denovos", main="Simulating De Novos in Well-Covered Regions", breaks=seq(min(three_or_more)-0.5, max(three_or_more)+0.5, 1), col="cyan", xaxt="n")
  axis(side=1, at = h3$mids, labels=seq(min(three_or_more),max(three_or_more),1))
  observed_triple = nrow(well_covered_regions[well_covered_regions$n_snp > 2,])
  abline(v=observed_triple, col="black", lty=3, lw=5)
  
}

plot_burden()

