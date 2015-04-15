# core routines for 'rupit' (regulatory mupit)

source("/Users//ps14/code//DDD/denovo_regs/R/regmut_null_model.R")

check_covered <- function(chr, pos, regions){
  
  # check if de novo is in regions
  
  chrom_regions = regions[regions$chr == chr, c("start", "stop")]
  match_pos = which(chrom_regions$start < pos & chrom_regions$stop > pos)
  if (length(match_pos) == 1){
    return(paste(chr, chrom_regions$start[match_pos], chrom_regions$stop[match_pos], sep = "."))
  } else {
    return(FALSE)
  }
}

check_contained <- function(chr, start, stop, regions){
  
  # check if region is CONTAINED in another region
  
  chrom_regions = regions[regions$chr == chr, c("start", "stop")]
  match_pos = which(chrom_regions$start <= start & chrom_regions$stop >= stop)
  if (length(match_pos) == 1){
    return(paste(chr, chrom_regions$start[match_pos], chrom_regions$stop[match_pos], sep = "."))
  } else {
    return(FALSE)
  }
}


filter_de_novos <- function(regions, de_novos){
  
  # filter de_novos to keep only:
    # 1. non-coding variants
    # 2. pp_dnm > 0.1 (denovogear prior)
    # 3. in specified regions (presumed to be well-covered)
  
  # adds region_id to keep track of which region de novo was located in
  
  noncoding = de_novos[de_novos$coding == "FALSE",]
  nc_filtered = noncoding[noncoding$pp_dnm > 0.1,]

  covered = unlist(mapply(check_covered, nc_filtered$chr, nc_filtered$pos, MoreArgs = list(regions)))
  
  de_novo_filtered = nc_filtered[covered != FALSE, ]
  de_novo_filtered$region_id = covered[covered != FALSE]
  return(de_novo_filtered)
}

generate_snp_null <- function(regions){
  
  # generates null model probability for snp at any point each region using trinucleotide mutation model
  # TODO: speed this code up
  
  return(unlist(lapply(regions$seq, p_sequence)))  # vector of p_nulls to be set as factor
}

generate_indel_null <- function(regions, indel2snp, global_mut_rate){
  
  # generates null model probability for indel in regions provided with global_mut_rate and indel2snp ratio over regions
  
  region_bp = sapply(as.character(regions$seq), nchar)
  total_bp = sum(region_bp)
  return(indel2snp*global_mut_rate*region_bp/total_bp)
}


count_denovos <- function(chr, start, stop, de_novos = TRUE){
  
  # count the number of snps and indels in each region that is passed through
  
  chr_match = de_novos[de_novos$chrom == chr,]
  pos_match = chr_match[(chr_match$pos > start) & (chr_match$pos < stop),]
  n_denovos = nrow(pos_match)
  if (n_denovos > 0){
    n_snp = nrow(pos_match[pos_match$var_type == "DENOVO-SNP",])
    n_indel = nrow(pos_match[pos_match$var_type == "DENOVO-INDEL",])
  } else {
    n_snp = 0
    n_indel = 0 }
  
  return(c(n_snp, n_indel))
}

assign_to_gene <- function(chr, start, stop, all_gencode){
  # reduce search to only genes on correct chromosome
  gencode_chrom = all_gencode[all_gencode$chr == paste("chr", chr, sep=""),] # paste to add chr
  
  # find closest start site accounting for strand
  distance_to_start = gencode_chrom$true_start - start
  distance_to_start[gencode_chrom$strand == "-"] = (-1)*distance_to_start[gencode_chrom$strand == "-"]
  
  distance_to_stop = gencode_chrom$true_stop - stop
  distance_to_stop[gencode_chrom$strand == "-"] = (-1)*distance_to_stop[gencode_chrom$strand == "-"]
  
  # check if start/stop are fully contained within an annotated gene. if so, assign this intronic region to that gene instead
  intronic = gencode_chrom[distance_to_start < 0 & distance_to_stop > 0,]
  if (nrow(intronic) == 1){
    closest_gene = intronic$gene
  } else {
    closest_start = which.min(sapply(distance_to_start, function(z) if (z > 0) z else Inf))
    closest_gene = gencode_chrom[closest_start, "gene"] }
  
  return(closest_gene)
}

regions_to_gencode <- function(regions, all_gencode){
  
  # return closest gene (by transcription start site)
  all_gencode <- all_gencode[all_gencode$gene_type == "protein_coding",]
  all_gencode$gene <- factor(all_gencode$gene)  # reset the factors (gene names)
  
  all_gencode$true_start <- all_gencode$start
  all_gencode$true_start[all_gencode$strand == "-"] <- all_gencode$stop[all_gencode$strand == "-"]
  
  all_gencode$true_stop <- all_gencode$stop
  all_gencode$true_stop[all_gencode$strand == "-"] <- all_gencode$start[all_gencode$strand == "-"]
  
  closest_genes = factor(mapply(assign_to_gene, regions$chr, regions$start, regions$stop, MoreArgs = list(all_gencode)))

  return(closest_genes)
}

get_gencode_sequence <- function(start, stop, region_id, regions){
  
  # adds $seq column to regions with sequence context taken from well_covered_regions (done by python script)
  
  region = regions[regions$region_id == region_id, c("start", "stop", "seq")]
  start_idx = start - region$start + 1
  stop_idx = stop - region$start
  return(substr(region$seq, start_idx, stop_idx))
  
  
}

expanded_regions <- function(regions, encode_dhs){
  
  # need to strip 'chr' off of dhs chromosome tag
  encode_dhs$dhs_chr = gsub("^chr","",encode_dhs$dhs_chr)
  
  contained = unlist(mapply(check_contained, encode_dhs$dhs_chr, encode_dhs$dhs_start, encode_dhs$dhs_end, MoreArgs = list(regions)))
  
  encode_dhs_contained = encode_dhs[contained != FALSE, ]
  encode_dhs_contained$region_id = contained[contained != FALSE]  # region id is ID for region it is CONTAINED IN - useful later on
  
  col_names = c("chr", "start", "stop", "closest_gene", "region_id", "seq")
  encode_dhs_contained$seq = unlist(mapply(get_gencode_sequence, encode_dhs_contained$dhs_start, encode_dhs_contained$dhs_end, encode_dhs_contained$region_id, MoreArgs = list(regions)))
  names(encode_dhs_contained)[names(encode_dhs_contained) %in% c("dhs_chr", "dhs_start", "dhs_end", "gene_name", "region_id", "seq")] = col_names
    
  regions = rbind(regions[,col_names], encode_dhs_contained[,col_names])
  return(regions)
}
