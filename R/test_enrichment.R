# test enrichment for de novo mutations in set of genomic regions. could be regulatory regions, or regulatory regions assigned
# to genes, for instance. grouping of regions should be specified by $region_id. could be gene name, or simple identifier

source("/Users/ps14/code/DDD/denovo_regs/R/rupit/R/rupit_core.R")

# load list of filtered de novos - run rupit_core.R function "filter_denovos" if not yet filtered
DENOVOS_FULL = "/Users/ps14/code/DDD/denovo_regs/R/rupit/data/de_novos.ddd_4k.noncoding_included.txt"
DENOVOS_FILTERED = "/Users/ps14/code/DDD/denovo_regs/R/rupit/data/de_novo_filtered.txt"
REGIONS_PATH = "/Users/ps14/code/DDD/denovo_regs/R/rupit/data/regions_annotated.txt"

# load data
de_novo_full <- read.table(DENOVOS_FULL, sep="\t", header=TRUE)
de_novo_filtered <- read.table(DENOVOS_FILTERED, sep="\t", header=TRUE)
well_covered_regions <- read.table(REGIONS_PATH, sep="\t", header=TRUE)

# check that REGIONS_PATH has $region_id - if not, assume each region to be tested independently and assign ids as chr.start.stop
if (!("region_id" %in% names(well_covered_regions))){
  well_covered_regions$region_id = paste(well_covered_regions$chr, well_covered_regions$start, well_covered_regions$stop, sep=".")
}

# check that null model priors have already been calculate
if (!("p_null" %in% names(well_covered_regions))){
  well_covered_regions$p_null = generate_prior(well_covered_regions)
}

# check that snps per region have been calculated. if not, update n_snp, n_indel
if (!("n_snp" %in% names(well_covered_regions))){
  well_covered_regions = update_counts(well_covered_regions, de_novo_filtered)
}

# set up a new dataframe 'regions' that groups regions by $region_id
reg = test_enrichment(well_covered_regions, by="region_id")

# set up a new dataframe 'genes' that groups regions associated with the same genes together
genes = test_enrichment(well_covered_regions, by="closest_gene")

# save enrichment analysis results
write.table(regions, file = "../data/regions_enrichment_pvals.txt", sep ="\t", col.names = TRUE)
write.table(genes, file = "../data/genes_enrichment_pvals.txt", sep ="\t", col.names = TRUE)



