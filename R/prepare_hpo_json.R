# script that prepares JSON files for region-wise de novo data
# these JSON files serve as input for HPO similarity calculation

# dependencies
library(rjson)

# load regions and filtered de novo regulatory variants (prepared by pre_process.R)
de_novo_filtered = read.table("../data/de_novo_filtered.txt", sep = "\t", header=TRUE)
well_covered_regions = read.table("../data/regions_annotated.txt", sep = "\t", header=TRUE)

# split into list with person_stable_id as key - this will be converted to a json
probands_by_region = split(de_novo_filtered$person_stable_id, de_novo_filtered$region_id)
region_json = toJSON(probands_by_region)
sink("../data/probands_by_region.json")
cat(region_json)
sink()
