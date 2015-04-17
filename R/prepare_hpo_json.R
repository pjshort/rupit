# script that prepares JSON files for region-wise de novo data
# these JSON files serve as input for HPO similarity calculation

# dependencies
library(jsonlite)

# load regions and filtered de novo regulatory variants (prepared by pre_process.R)
de_novo_filtered = read.table("../data/de_novo_filtered.txt", sep = "\t", header=TRUE)
well_covered_regions = read.table("../data/regions_annotated.txt", sep = "\t", header=TRUE)

# remove diagnosed probands - comment out if not needed
# TODO - refactor to make this better...
# diagnosed_probands <- read.table("../data/ddd_likely_diagnosed.txt", sep="\t", header=TRUE)
# de_novo_filtered <- de_novo_filtered[!(de_novo_filtered$person_stable_id %in% diagnosed_probands$person_id),]

# split into list with person_stable_id as key - this will be converted to a json
split_to_json <- function(de_novos, by){
  probands_by_factor = split(de_novos$person_stable_id, de_novos[, by])
  probands_by_factor = Filter(function(x) length(x) > 0, probands_by_factor)
  region_json = toJSON(probands_by_factor, pretty = TRUE)
  save_name = sprintf("../data/probands_by_%s.json", by)
  sink(save_name)
  cat(region_json)
  sink()
}

# split probands by region_id
split_to_json(de_novo_filtered, by = "region_id")

# split probands by closest_gene
de_novo_filtered = merge(de_novo_filtered, well_covered_regions[,c("region_id", "closest_gene")], by = "region_id")
split_to_json(de_novo_filtered, by = "closest_gene")

# split probands by regulatory region annotation summary (heart, conserved, enhancer, combinations...)
#de_novo_filtered = merge(de_novo_filtered, well_covered_regions[,c("region_id", "annotsummary")], by = "region_id")
#split_to_json(de_novo_filtered, by = "annotsummary")






