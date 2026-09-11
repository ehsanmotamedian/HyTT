# 5' UTR ensemble free energy, MEA structure/energy, and mean base-pair
# probability via ViennaRNA RNAfold -p --MEA
# Requires: seqinr, openxlsx, and a working ViennaRNA installation with
# RNAfold available on PATH.

# install.packages(c("seqinr", "openxlsx"))

library(seqinr)
library(openxlsx)

data_dir   <- "data"
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

orfs <- readLines(file.path(data_dir, "gene_list.txt"))

rnafold_cmd <- Sys.which("RNAfold")
if (rnafold_cmd == "") rnafold_cmd <- Sys.which("RNAfold.exe")
if (rnafold_cmd == "") {
  stop("RNAfold not found on PATH. Install ViennaRNA and ensure RNAfold is accessible.")
}

fasta_file <- file.path(data_dir, "SGD_all_ORFs_5prime_UTRs.fsa")
sequences_list <- read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE)

sequences <- data.frame(gene_id = character(), seq = character(), stringsAsFactors = FALSE)
for (orf in orfs) {
  matching_seq <- grep(orf, names(sequences_list), value = TRUE)
  if (length(matching_seq) > 0) {
    seq_str <- as.character(sequences_list[[matching_seq[1]]])
    seq_str <- gsub("T", "U", seq_str)
    sequences <- rbind(sequences, data.frame(gene_id = orf, seq = seq_str, stringsAsFactors = FALSE))
  } else {
    sequences <- rbind(sequences, data.frame(gene_id = orf, seq = "", stringsAsFactors = FALSE))
  }
}

results <- data.frame(
  ORF = orfs,
  Ensemble_Free_Energy = NA,
  Centroid_Structure = NA,
  MEA_Structure = NA,
  MEA_Energy = NA,
  Mean_Base_Pair_Prob = NA,
  stringsAsFactors = FALSE
)

parse_rnafold_output <- function(output) {
  ensemble_energy <- NA
  centroid_struct <- NA
  mea_struct <- NA
  mea_energy <- NA

  ensemble_line <- grep("\\[", output, value = TRUE)
  if (length(ensemble_line) > 0) {
    m <- regmatches(ensemble_line[1], regexec("\\[\\s*([-0-9.]+)\\]", ensemble_line[1]))
    if (length(m[[1]]) > 1) ensemble_energy <- as.numeric(m[[1]][2])
  }

  centroid_line <- grep("d=", output, value = TRUE)
  if (length(centroid_line) > 0) {
    centroid_struct <- gsub("^\\s*([.()]+)\\s*.*", "\\1", centroid_line[1])
  }

  mea_line <- grep("MEA=", output, value = TRUE)
  if (length(mea_line) > 0) {
    mea_struct <- gsub("^\\s*([.()]+)\\s*.*", "\\1", mea_line[1])
    m <- regmatches(mea_line[1], regexec("MEA=\\s*([0-9.]+)", mea_line[1]))
    if (length(m[[1]]) > 1) mea_energy <- as.numeric(m[[1]][2])
  }

  return(list(ensemble_energy = ensemble_energy, centroid_struct = centroid_struct,
              mea_struct = mea_struct, mea_energy = mea_energy))
}

for (i in seq_along(orfs)) {
  seq <- sequences$seq[sequences$gene_id == orfs[i]]
  if (length(seq) == 0 || nchar(seq) == 0) {
    cat("No sequence found for", orfs[i], "\n")
    next
  }
  # Cap sequence length at 500 nt to keep the partition-function (-p) calculation tractable
  if (nchar(seq) > 500) {
    seq <- substr(seq, 1, 500)
  }

  temp_dir <- tempfile(pattern = "rnafold-")
  dir.create(temp_dir)
  old_wd <- getwd()
  setwd(temp_dir)

  temp_file <- "input.fasta"
  writeLines(c(paste(">", orfs[i]), seq), temp_file)

  cmd <- paste(shQuote(rnafold_cmd), "-p --MEA", temp_file)
  output <- system(cmd, intern = TRUE)

  setwd(old_wd)

  if (length(output) == 0 || grepl("ERROR", output[1])) {
    cat("RNAfold error for", orfs[i], "\n")
    unlink(temp_dir, recursive = TRUE)
    next
  }

  parsed <- parse_rnafold_output(output)
  results$Ensemble_Free_Energy[i] <- parsed$ensemble_energy
  results$Centroid_Structure[i] <- parsed$centroid_struct
  results$MEA_Structure[i] <- parsed$mea_struct
  results$MEA_Energy[i] <- parsed$mea_energy

  dp_file <- file.path(temp_dir, paste0(orfs[i], "_dp.ps"))
  if (file.exists(dp_file)) {
    dp_lines <- readLines(dp_file)
    ubox_lines <- grep("ubox$", dp_lines, value = TRUE)
    sum_p <- 0
    if (length(ubox_lines) > 0) {
      for (line in ubox_lines) {
        line <- trimws(line)
        parts <- strsplit(line, "\\s+")[[1]]
        if (length(parts) >= 4 && parts[4] == "ubox") {
          sqp <- as.numeric(parts[3])
          if (!is.na(sqp)) {
            sum_p <- sum_p + sqp^2
          }
        }
      }
    }
    seq_len <- nchar(seq)
    results$Mean_Base_Pair_Prob[i] <- if (seq_len > 0) 2 * sum_p / seq_len else NA

    unlink(dp_file)
    unlink(file.path(temp_dir, paste0(orfs[i], "_ss.ps")))
  } else {
    cat("Dot-plot file not generated for", orfs[i], "\n")
  }

  unlink(temp_dir, recursive = TRUE)
}

write.xlsx(results, file.path(output_dir, "ensemble_free_energy.xlsx"))
