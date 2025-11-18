
# Load relevant packages

library(jsonlite)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(openxlsx)


# 1. Load JSON FILE


json_file <- "C:/Users/sulsowe/Documents/WORKINGS/results_Batch11.json"
data <- fromJSON(json_file, simplifyVector = FALSE)

samples <- data$samples


# 2. Define virulence gene groups 


virulence_groups <- list(
  clf = c("clfA","clfB"),
  fnb = c("fnbA","fnbB"),
  spa = "spa",
  cna = "cna",
  ebpS = "ebpS",
  coa_vwb = c("coa","vwb"),
  sak = "sak",
  hysA = "hysA",
  lip = "lip",
  nuc = "nuc",
  hla_hly = c("hla","hly"),
  hlb = "hlb",
  hld = "hld",
  hlg = "hlg",
  pvl = c("lukS","lukF"),
  tst = "tst",
  sea_see = c("sea","seb","sec","sed","see","seg","seh","sei","sej"),
  eta_etb = c("eta","etb"),
  cap = c("cap","capA","capB","capC","capD","capE","capF","capG","capH","capI","capJ","capK"),
  chp = "chp",
  isd = c("isdA","isdB","isdC"),
  sfa_sfb = c("sfa","sfb"),
  psm = c("psmA","psmB","psmα","psmβ"),
  agr = c("agrA","agrB","agrC","agrD"),
  sarA = "sarA",
  sae = c("saeR","saeS"),
  ica = c("icaA","icaB","icaC","icaD"),
  bap = "bap",
  aur = "aur",
  lta = "lipoteichoic acid",
  peptidoglycan = "peptidoglycan"
)


# Helper function: detect virulence gene groups


detect_virulence <- function(annotation_list, groups) {
  
  if (is.null(annotation_list) || length(annotation_list) == 0) {
    return(setNames(rep(FALSE, length(groups)), names(groups)))
  }
  
  annotation_text <- tolower(
    paste(
      unlist(lapply(annotation_list, function(x) {
        paste(x$gene %||% "", x$product %||% "")
      })),
      collapse = " "
    )
  )
  
  sapply(groups, function(genes) {
    any(stringr::str_detect(annotation_text, tolower(genes)))
  })
}

# AMR Extraction

extract_amr <- function(smp) {
  amr_hits <- smp$results$antimicrobial_resistance$detected_variants
  
  if (is.null(amr_hits) || length(amr_hits) == 0) {
    return(character(0))
  }
  
  unique(vapply(amr_hits, function(x) x$gene, character(1)))
}


# Extract data for each sample

result_list <- lapply(samples, function(smp) {
  
  sample_id <- smp$alias
  
  st <- smp$results$sequence_typing$sequence_type %||% NA
  
  amr_genes <- extract_amr(smp)
  
  annotations <- smp$results$assembly$annotations
  
  virulence <- detect_virulence(annotations, virulence_groups)
  
  list(
    sample = sample_id,
    ST = st,
    AMR = amr_genes,
    virulence = virulence
  )
})



# Build AMR yes/no matrix


all_amr_genes <- sort(unique(unlist(lapply(result_list, `[[`, "AMR"))))

amr_matrix <- do.call(rbind, lapply(result_list, function(x) {
  as.data.frame(t(all_amr_genes %in% x$AMR))
}))

colnames(amr_matrix) <- all_amr_genes


# Build virulence yes/no matrix


vir_matrix <- do.call(rbind, lapply(result_list, function(x) {
  as.data.frame(t(x$virulence))
}))


# Convert T/F to Yes/No

logical_to_yesno <- function(df) {
  df[] <- lapply(df, function(x) ifelse(x, "Yes", "No"))
  df
}

amr_matrix <- logical_to_yesno(amr_matrix)
vir_matrix <- logical_to_yesno(vir_matrix)

# Combine all into final table

final_df <- data.frame(
  Sample = sapply(result_list, `[[`, "sample"),
  ST = sapply(result_list, `[[`, "ST"),
  amr_matrix,
  vir_matrix,
  check.names = FALSE,
  stringsAsFactors = FALSE
)



# Save to Excel (Optional)


write.xlsx(final_df, "Batch11_json_results.xlsx", overwrite = TRUE)

