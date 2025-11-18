# Load packages

library(dplyr)
library(stringr)

# 1. Folder containing all GFF files

gff_folder <- "C:/Users/sulsowe/Documents/WORKINGS/gff_files_Batch12"

gff_files <- list.files(gff_folder, pattern = "\\.gff$", full.names = TRUE)


# 2. Virulence gene list (your 30 virulence groups)

virulence_list <- list(
  clfA_clfB = c("clfA", "clfB"),
  fnbA_fnbB = c("fnbA", "fnbB"),
  spa = "spa",
  cna = "cna",
  ebpS = "ebpS",
  coa_vwb = c("coa", "vwb"),
  sak = "sak",
  hysA = "hysA",
  lip = "lip",
  nuc = "nuc",
  hla_hly = c("hla", "hly"),
  hlb = "hlb",
  hld = "hld",
  hlg = "hlg",
  lukPV = c("lukS", "lukF"),
  tst = "tst",
  sea_seg_sei = c("sea", "seb", "sec", "sed", "see", "seg", "seh", "sei"),
  eta_etb = c("eta", "etb"),
  cap_operon = c("capA","capB","capC","capD","capE","capF","capG","capH","capI","capJ","capK","capL","capM","capN","capO"),
  chp = "chp",
  isdABC = c("isdA","isdB","isdC"),
  sfa_sfb = c("sfa","sfb"),
  psm_alpha_beta = c("psmA", "psmB"),
  agrABCD = c("agrA","agrB","agrC","agrD"),
  sarA = "sarA",
  saeRS = c("saeR","saeS"),
  icaABCD = c("icaA","icaB","icaC","icaD"),
  bap = "bap",
  aur = "aur",
  lipoteichoic_acid = "lipoteichoic",
  peptidoglycan = "peptidoglycan"
)


# 3. Function to read GFF and search virulence genes


scan_gff <- function(gff_path) {
  
  # Read only non-comment lines (skip ## lines)
  gff <- read.delim(
    gff_path,
    header = FALSE,
    comment.char = "#",
    sep = "\t",
    quote = "",
    stringsAsFactors = FALSE,
    fill = TRUE
  )
  
  # Ensure we have columns
  if (ncol(gff) < 9) return(NULL)
  
  colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")
  
  # Extract gene and product from the attributes column
  attributes <- tolower(gff$attributes)
  
  gene_matches <- str_extract(attributes, "gene=[^;]+")
  product_matches <- str_extract(attributes, "product=[^;]+")
  
  gene_vals <- tolower(str_replace(gene_matches, "gene=", ""))
  prod_vals <- tolower(str_replace(product_matches, "product=", ""))
  
  combined <- unique(c(gene_vals, prod_vals))
  combined <- combined[!is.na(combined)]
  
  # Test each virulence group
  detected <- sapply(virulence_list, function(glist) {
    any(sapply(glist, function(g) any(str_detect(combined, paste0("\\b", g, "\\b")))))
  })
  
  return(detected)
}


# 4. Apply to every GFF file

vir_data <- lapply(gff_files, scan_gff)
names(vir_data) <- basename(gff_files)

vir_matrix <- do.call(rbind, vir_data)

# Convert TRUE/FALSE → Yes/No
vir_matrix <- as.data.frame(apply(vir_matrix, 2, function(x) ifelse(x, "Yes", "No")))

# Add sample names
vir_matrix <- cbind(Sample = rownames(vir_matrix), vir_matrix)
rownames(vir_matrix) <- NULL

