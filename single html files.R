

# Relevant packages

library(rvest)
library(dplyr)
library(purrr)
library(stringr)
library(openxlsx)

# 1. Set working directory containing all single html barcode files

folder <- "C:/Users/sulsowe/Documents/WORKINGS/Batch14_single_html_files"   # INPUT FOLDER
files <- list.files(folder, pattern = "\\.html?$", full.names = TRUE)

# VIRULENCE GENE GROUPS TO SCREEN

virulence_groups <- list(
  clfA_clfB            = c("clfA","clfB"),
  fnbA_fnbB            = c("fnbA","fnbB"),
  spa                  = c("spa"),
  cna                  = c("cna"),
  ebpS                 = c("ebpS"),
  coa_vwb              = c("coa","vwb"),
  sak                  = c("sak"),
  hysA                 = c("hysA"),
  lip                  = c("lip"),
  nuc                  = c("nuc"),
  hla_hly              = c("hla","hly"),
  hlb                  = c("hlb"),
  hld                  = c("hld"),
  hlg_lukSPV_lukFPV    = c("hlg", "lukS-PV", "lukF-PV"),
  tst                  = c("tst"),
  sea_to_sei           = c("sea","seb","sec","sed","see","seg","seh","sei"),
  eta_etb              = c("eta","etb"),
  cap_operon           = c("cap"),
  chp                  = c("chp"),
  isdA_B_C             = c("isdA","isdB","isdC"),
  sfa_sfb              = c("sfa","sfb"),
  psm_alpha_beta       = c("psmα","psmβ","psmA","psmB"),
  agr_operon           = c("agrA","agrB","agrC","agrD"),
  sarA                 = c("sarA"),
  saeR_S               = c("saeR","saeS"),
  icaA_D_operon        = c("icaA","icaB","icaC","icaD"),
  bap                  = c("bap"),
  aur                  = c("aur"),
  lipoteichoic_acid    = c("lipoteichoic acid"),
  peptidoglycan        = c("peptidoglycan")
)

# MLST extraction using XPath

extract_mlst <- function(html){
  mlst_node <- html %>% 
    html_element(xpath = "//h2[contains(., 'MLST')]//following::table[1]")
  
  if (is.na(mlst_node)) return(NULL)
  
  mlst_table <- mlst_node %>% html_table()
  
  # convert all MLST columns to character to avoid mixing integer/character
  mlst_table <- mlst_table %>% mutate(across(everything(), as.character))
  
  return(mlst_table)
}


# AMR extraction

extract_amr <- function(html){
  nodes <- html %>% 
    html_elements(xpath = "//h2[contains(., 'Antimicrobial')]/following::button//b")
  html_text(nodes)
}


# VIRULENCE extraction from ANNOTATIONS

extract_annotations <- function(html){
  ann <- html %>% 
    html_element(xpath = "//h2[contains(., 'Annotation')]//following::table[1]")
  
  if (is.na(ann)) return(character())
  
  ann_table <- ann %>% html_table()
  
  # Combine all annotation text into one searchable string
  annotation_text <- paste(ann_table$Product, ann_table$Gene, collapse="; ", sep=" ")
  annotation_text <- tolower(annotation_text)
  
  return(annotation_text)
}


# Screen for virulence genes in annotation text

detect_virulence <- function(annotation_text){
  results <- sapply(virulence_groups, function(glist){
    
    found <- any(sapply(glist, function(g)
      grepl(tolower(g), annotation_text, fixed = TRUE)
    ))
    
    ifelse(found, "yes", "no")
  })
  
  return(results)
}


# Process each file

all_results <- map(files, function(f){
  
  html <- read_html(f)
  sample <- tools::file_path_sans_ext(basename(f))
  
  mlst <- extract_mlst(html)
  amr  <- extract_amr(html)
  ann  <- extract_annotations(html)
  vir  <- detect_virulence(ann)
  
  if (!is.null(mlst)) {
    mlst$Sample <- sample
  }
  
  list(sample = sample, mlst = mlst, amr = amr, vir = vir)
})


# AMR yes/no matrix

all_amr_genes <- unique(unlist(map(all_results, "amr")))

amr_matrix <- map_dfr(all_results, function(x){
  tibble(
    Sample = x$sample,
    !!!setNames(lapply(all_amr_genes, function(g) 
      ifelse(g %in% x$amr, "yes", "no")), all_amr_genes)
  )
})


# Virulence yes/no matrix

vir_matrix <- map_dfr(all_results, function(x){
  tibble(
    Sample = x$sample,
    !!!as.list(x$vir)
  )
})


# Combine everything

mlst_table <- bind_rows(map(all_results, "mlst"))

final <- mlst_table %>%
  left_join(amr_matrix, by = "Sample") %>%
  left_join(vir_matrix, by = "Sample")

# Export to excel (optional)

write.xlsx(final, "batch14_epi2me_results.xlsx")

