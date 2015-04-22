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
iterations = 100

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

simulate_de_novos_by_gene <- function(regions, snp_total, iterations){
  
  # regions are aggregated by the $closest_gene then de novos are simulated in the exact same manner
  
  p_null_gene = aggregate(well_covered_regions$p_snp_null, by = list(closest_gene = well_covered_regions$closest_gene), FUN = sum)
  genes = do.call(cbind.data.frame, p_null_gene)
  names(genes)[names(genes) == "x"] = "p_snp_null"
  
  sim_output = simulate_de_novos(genes, snp_total, iterations)
  return(sim_output)
}

plot_burden <- function(sim_results, observed_hits){
  
  # plots histograms of the number of genomic regions with 1 de novo, 2+ de novos
  # takes results of simulation (large matrix) and observed hits (length 2 vector corresp to # regions with ==1 and >1)
    
  single_hits = colSums(sim_results == 1)
  two_or_more = colSums(sim_results > 1)
  
  # single de novos
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  h1 = hist(single_hits, xlab="# of regions with 1 denovo", main = "Single De Novo", breaks=seq(min(single_hits) - 0.5, max(single_hits)+0.5, 1), col="cyan", xaxt="n")
  axis(side=1, at = h1$mids, labels=seq(min(single_hits),max(single_hits),1))
  abline(v=observed_hits[1], col="black", lty=3, lw=5)
  
  # double de novos
  par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
  h2 = hist(two_or_more, xlab="# of regions with 2+ denovos", main = "Recurrent De Novos", breaks=seq(min(two_or_more)-0.5, max(two_or_more)+0.5, 1), col="cyan", xaxt="n")
  axis(side=1, at = h2$mids, labels=seq(min(two_or_more),max(two_or_more),1))
  observed_double = nrow(well_covered_regions[well_covered_regions$n_snp == 2,])
  abline(v=observed_hits[2], col="black", lty=3, lw=5)
  
}

simulate_gene <- function(gene_name, sim_results, plot = FALSE, observed_hits = NULL){
  
  # return results of simulation for a single gene. plots # observed mutations as histogram plot = TRUE
  
  region_match = which(well_covered_regions$closest_gene == gene_name)
  gene_results = sim_results[region_match,]
  if (length(region_match) > 1){
    counts = colSums(gene_results)
  } else {
    counts = gene_results
  }
  
  if (plot == TRUE){
    par(mar=c(5,4,5,1) + 0.5)   # extra large bottom margin
    h = hist(counts, xlab="# of hits in gene-associated regulatory regions", main=paste0("Simulating Recurrent De Novos in ", gene_name), breaks=seq(min(counts)-0.5, max(max(counts), observed_hits)+0.5, 1), col="cyan", xaxt="n")
    if (!is.null(observed_hits)){
      abline(v=observed_hits, col="black", lty=3, lw=5)
      axis(side=1, at = seq(min(h$mids), max(max(h$mids), max(observed_hits))), labels=seq(min(counts),max(max(counts), observed_hits),1))
    } else {
      axis(side=1, at = h$mids, labels=seq(min(counts),max(counts),1))
    }
  }

  return(gene_results)
}

gene_sim_likelihood <- function(gene_name, sim_results, observed_hits){
  
  # calculates the probability that recurrent de novos >2 number of observed hits were observed in simulation
  # this is a proxy for likelihood of the observed event under the null model
  
  region_match = which(well_covered_regions$closest_gene == gene_name)
  gene_results = sim_results[region_match,]
  if (length(region_match) > 1){
    counts = colSums(gene_results)
  } else {
    counts = gene_results
  }

  trials = length(counts)
  success = sum(counts >= observed_hits)
  return(success/trials)
  
}


