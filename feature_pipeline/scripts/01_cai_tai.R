# CAI and tAI calculation for S. cerevisiae ORFs
# Requires: biomaRt, seqinr, tAI, jsonlite, openxlsx

# install.packages(c("seqinr", "tAI", "jsonlite", "openxlsx"))
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("biomaRt")

library(biomaRt)
library(seqinr)
library(tAI)
library(jsonlite)
library(openxlsx)

data_dir   <- "feature_pipeline/data"
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

# ---- tRNA gene copy numbers for S. cerevisiae (GtRNAdb; Chan & Lowe, 2016) ----
# Vector length 64 (standard codon order), used directly by tAI::get.ws.
trna <- rep(0, 64)
trna[3]  <- 14 # AAG: Lys-CTT
trna[1]  <- 1  # AAA: Lys-TTT
trna[2]  <- 10 # AAC: Asn-GTT
trna[8]  <- 11 # ACT: Thr-AGT
trna[7]  <- 1  # ACG: Thr-CGT
trna[5]  <- 1  # ACA: Thr-TGT
trna[24] <- 10 # CCT: Pro-AGG
trna[23] <- 1  # CCG: Pro-CGG
trna[21] <- 1  # CCA: Pro-TGG
trna[9]  <- 6  # AGA: Arg-TCT
trna[11] <- 4  # AGG: Arg-CCT
trna[28] <- 5  # CGT: Arg-ACG
trna[27] <- 1  # CGG: Arg-CCG
trna[16] <- 13 # ATT: Ile-AAT
trna[13] <- 1  # ATA: Ile-TAT
trna[15] <- 13 # ATG: Met-CAT
trna[19] <- 4  # CAG: Gln-CTG
trna[17] <- 9  # CAA: Gln-TTG
trna[18] <- 4  # CAC: His-GTG
trna[37] <- 5  # GCA: Ala-UGC
trna[40] <- 11 # GCT: Ala-AGC
trna[42] <- 16 # GGC: Gly-GCC
trna[41] <- 5  # GGA: Gly-TCC
trna[48] <- 12 # GTT: Val-AAC
trna[47] <- 2  # GTG: Val-CAC
trna[45] <- 2  # GTA: Val-TAC
trna[59] <- 3  # TGG: Trp-CCA
trna[49] <- 0  # TAA: Stop
trna[50] <- 8  # TAC: Tyr-GTA
trna[51] <- 0  # TAG: Stop
trna[56] <- 11 # TCT: Ser-AGA
trna[55] <- 3  # TCG: Ser-CGA
trna[53] <- 1  # TCA: Ser-TGA
trna[10] <- 5  # AGC: Ser-GCT
trna[57] <- 0  # TGA: Stop
trna[58] <- 4  # TGC: Cys-GCA
trna[32] <- 11 # CTT: Leu-AAG
trna[63] <- 10 # TTG: Leu-CAA
trna[31] <- 1  # CTG: Leu-CAG
trna[61] <- 1  # TTA: Leu-TAA
trna[29] <- 1  # CTA: Leu-TAG
trna[62] <- 10 # TTC: Phe-GAA
trna[34] <- 16 # GAC: Asp-GTC
trna[35] <- 14 # GAG: Glu-CTC
trna[33] <- 1  # GAA: Glu-TTC

# ---- Retrieve CDS from SGD (biomaRt is unreliable for this endpoint) ----
get_cds <- function(orf_name) {
  url <- paste0("https://www.yeastgenome.org/backend/locus/", orf_name, "/sequence_details")
  response <- fromJSON(url, simplifyVector = FALSE)
  coding_residues <- NULL
  if ("coding_dna" %in% names(response)) {
    for (seq_obj in response$coding_dna) {
      if (seq_obj$strain$display_name == "S288C") {
        coding_residues <- seq_obj$residues
        break
      }
    }
  }
  if (is.null(coding_residues)) {
    stop("Coding sequence not found in SGD response for S288C.")
  }
  return(s2c(coding_residues))
}

calculate_cai <- function(cds_seq) {
  data(caitab) # seqinr's built-in codon usage weight table
  w <- caitab$sc # S. cerevisiae weights
  cai(cds_seq, w = w, numcode = 1)
}

calculate_tai <- function(cds_seq) {
  ws <- get.ws(tRNA = trna, sking = 0) # sking = 0 for eukaryotes
  codon_freq <- uco(cds_seq, index = "f", as.data.frame = FALSE)
  # Exclude Met and stop codons (dos Reis et al., 2004)
  stop_met_pos <- c(15, 49, 51, 57) # ATG, TAA, TAG, TGA
  codon_freq[stop_met_pos] <- 0
  codon_freq <- codon_freq / sum(codon_freq)
  codon_freq <- codon_freq[-stop_met_pos]
  codon_freq <- matrix(codon_freq, nrow = 1)
  get.tai(codon_freq, ws)
}

# ---- Gene list ----
orf_names <- readLines(file.path(data_dir, "gene_list.txt"))

# ---- Compute CAI and tAI for each gene ----
results <- data.frame(ORF = character(), CAI = numeric(), tAI = numeric(), Error = character(), stringsAsFactors = FALSE)

for (orf_name in orf_names) {
  tryCatch({
    cds_seq <- get_cds(orf_name)
    cai_val <- calculate_cai(cds_seq)
    tai_val <- calculate_tai(cds_seq)
    results <- rbind(results, data.frame(ORF = orf_name, CAI = cai_val, tAI = tai_val, Error = NA))
    cat(orf_name, "- CAI:", cai_val, "- tAI:", tai_val, "\n")
  }, error = function(e) {
    results <<- rbind(results, data.frame(ORF = orf_name, CAI = NA, tAI = NA, Error = e$message))
    cat(orf_name, "- Error:", e$message, "\n")
  })
}

write.xlsx(results, file.path(output_dir, "cai_tai_results.xlsx"), rowNames = FALSE)
