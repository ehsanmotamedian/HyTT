# 5' and 3' UTR minimum free energy (MFE) via ViennaRNA RNAfold
# Requires: seqinr, stringr, openxlsx, and a working ViennaRNA installation
# (https://www.tbi.univie.ac.at/RNA/) with RNAfold available on PATH.

# install.packages(c("seqinr", "stringr", "openxlsx"))

library(seqinr)
library(stringr)
library(openxlsx)

data_dir   <- "data"
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

genes <- readLines(file.path(data_dir, "gene_list.txt"))

utr5_path <- file.path(data_dir, "SGD_all_ORFs_5prime_UTRs.fsa")
utr3_path <- file.path(data_dir, "SGD_all_ORFs_3prime_UTRs.fsa")

# Locate the RNAfold executable (works whether it's called RNAfold or RNAfold.exe)
rnafold_cmd <- Sys.which("RNAfold")
if (rnafold_cmd == "") rnafold_cmd <- Sys.which("RNAfold.exe")
if (rnafold_cmd == "") {
  stop("RNAfold not found on PATH. Install ViennaRNA and ensure RNAfold is accessible.")
}

extract_sequences <- function(file_path, genes) {
  sequences <- read.fasta(file_path, seqtype = "DNA", as.string = TRUE)
  seqs <- vector("list", length(genes))
  names(seqs) <- genes
  for (i in seq_along(genes)) {
    gene <- genes[i]
    matching <- grep(toupper(gene), toupper(names(sequences)), fixed = TRUE)
    if (length(matching) > 0) {
      seqs[[gene]] <- toupper(sequences[[matching[1]]])
      seqs[[gene]] <- gsub("T", "U", seqs[[gene]]) # DNA -> RNA
    } else {
      seqs[[gene]] <- NA
    }
  }
  return(seqs)
}

utr5_seqs <- extract_sequences(utr5_path, genes)
utr3_seqs <- extract_sequences(utr3_path, genes)

calculate_mfe <- function(sequence) {
  if (is.na(sequence) || nchar(sequence) == 0) {
    return(NA)
  }
  temp_input <- tempfile(fileext = ".fa")
  writeLines(c(">seq", sequence), temp_input)
  cmd <- paste0(shQuote(rnafold_cmd), " --noPS --infile=", shQuote(temp_input))
  result <- system(cmd, intern = TRUE)
  mfe_line <- tail(result, 1)
  mfe <- as.numeric(gsub(".*\\(\\s*([-+]?[0-9]*\\.?[0-9]+)\\s*\\).*", "\\1", mfe_line))
  unlink(temp_input)
  return(mfe)
}

utr5_mfe <- sapply(utr5_seqs, calculate_mfe)
utr3_mfe <- sapply(utr3_seqs, calculate_mfe)

df <- data.frame(Gene = genes, UTR5_MFE = utr5_mfe, UTR3_MFE = utr3_mfe)
write.xlsx(df, file.path(output_dir, "mfe_results.xlsx"), rowNames = FALSE)
