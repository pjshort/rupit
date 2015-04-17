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

# TODO - determine a better way to toggle diagnosed/undiagnosed on and off
#diagnosed_probands <- read.table("../data/ddd_likely_diagnosed.txt", sep="\t", header=TRUE)
#de_novo_full <- de_novo_full[!(de_novo_full$person_stable_id %in% diagnosed_probands$person_id),]

# check that REGIONS_PATH has $region_id - if not, assume each region to be tested independently and assign ids as chr.start.stop
if (any(!(c("region_id", "p_snp_null", "p_indel_null", "n_snp", "n_indel") %in% names(well_covered_regions)))){
  print("Need to run pre_process.R on regions, de novos first!")
}

# set up a new dataframe 'regions' that groups regions by $region_id
regions = test_enrichment(well_covered_regions, by="region_id")
regions = merge(regions, well_covered_regions[,c("region_id", "closest_gene")], by = "region_id")


# set up a new dataframe 'genes' that groups regions associated with the same genes together
genes = test_enrichment(well_covered_regions, by="closest_gene")
names(genes)[names(genes) == "region_id"] = "closest_gene"

# save enrichment analysis results
write.table(regions, file = "../data/regions_enrichment_pvals.txt", sep ="\t", col.names = TRUE)
write.table(genes, file = "../data/genes_enrichment_pvals.txt", sep ="\t", col.names = TRUE)



