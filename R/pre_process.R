# pre-process data to greatly speed up analysis

# if starting from SCRATCH with just regions, de novos, this is the first thing that should be run.
# it may take a while, but will greatly speed up later analysis

# requirements:
  # 1. de novo mutations file (TSV) with chr, pos, pp_dnmm, coding column names required
  # 2. genomic regions with chr, start, stop, region_id (if planning to test for enrichment later - i.e. the closest gene)

# output files:
  # de novos filtered down to include only those that lie in regions provided, with pp_dnm > 0.1
  # regions annotated with n_snp (from de novo file), p_null (prior prob of mutation from null model)
  # if ADD_GENCODE == TRUE, $closest_gene will be added to regions file

# filenames:
  # de_novo_filtered.txt
  # regions_annotated.txt

source("/Users/ps14/code/DDD/denovo_regs/R/rupit/R/rupit_core.R")

# files will be saved to rupit/data unless location is changed/specified:
OUT_PATH = "/Users/ps14/code/DDD/denovo_regs/R/rupit/data"

ADD_GENCODE = TRUE
GENCODE_PATH = "/Users/ps14/code/DDD/denovo_regs/R/rupit/data/all_gencode_genes_v19_+strand.txt"

DENOVOS_PATH = "/Users/ps14/code/DDD/denovo_regs/R/rupit/data/de_novos.ddd_4k.noncoding_included.txt"
REGIONS_PATH = "/Users/ps14/code/DDD/denovo_regs/R/rupit/data/DDD_TRRs.annotated.highcov.sequence.txt"

de_novo_full <- read.table(DENOVOS_PATH, sep="\t", header=TRUE)
well_covered_regions <- read.table(REGIONS_PATH, sep="\t", header=TRUE)

# call filter functions from rupit_core
de_novo_filtered = filter_de_novos(well_covered_regions, de_novo_full)
well_covered_regions = update_counts(well_covered_regions)
well_covered_regions$region_id = paste(well_covered_regions$chr, well_covered_regions$start, well_covered_regions$stop, sep = ".")

well_covered_regions$p_null = generate_prior(well_covered_regions)

if (ADD_GENCODE == TRUE){  
  all_gencode <- read.table(GENCODE_PATH, sep="\t", header=TRUE)
  well_covered_regions = regions_to_gencode(well_covered_regions, all_gencode)
}

write.table(de_novo_filtered, file = paste0(OUT_PATH,"/de_novo_filtered.txt"), sep = "\t")
write.table(well_covered_regions, file = paste0(OUT_PATH,"/regions_annotated.txt"), sep = "\t")
            
