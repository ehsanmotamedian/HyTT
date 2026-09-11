# 5' and 3' UTR lengths for S. cerevisiae genes
# Requires: seqinr, stringr, openxlsx

# install.packages(c("seqinr", "stringr", "openxlsx"))

library(seqinr)
library(stringr)
library(openxlsx)

data_dir   <- "data"
output_dir <- "output"
dir.create(output_dir, showWarnings = FALSE)

genes <- readLines(file.path(data_dir, "gene_list.txt"))

# Extract UTR length per gene; averages across isoforms when a gene has more than one
extract_utr_lengths_for_genes <- function(file_path, genes) {
  sequences <- read.fasta(file_path, seqtype = "DNA", as.string = FALSE)
  lengths <- vector("numeric", length(genes))
  for (i in seq_along(genes)) {
    gene <- genes[i]
    matching <- grep(gene, names(sequences), fixed = TRUE)
    if (length(matching) > 0) {
      lens <- sapply(sequences[matching], length)
      lengths[i] <- mean(lens, na.rm = TRUE)
    } else {
      lengths[i] <- NA
    }
  }
  return(lengths)
}

utr5_lengths <- extract_utr_lengths_for_genes(file.path(data_dir, "SGD_all_ORFs_5prime_UTRs.fsa"), genes)
utr3_lengths <- extract_utr_lengths_for_genes(file.path(data_dir, "SGD_all_ORFs_3prime_UTRs.fsa"), genes)

df <- data.frame(Gene = genes, UTR5_Length = utr5_lengths, UTR3_Length = utr3_lengths)
write.xlsx(df, file.path(output_dir, "utr_lengths.xlsx"), rowNames = FALSE)
