# Feature Engineering Pipeline

Eighteen scripts compute 115 sequence-derived, transcriptomic, and biophysical
features for 1710 *S. cerevisiae* genes. Run in the numeric order below.

## One-time setup

1. Set a valid email address for NCBI Entrez queries. Each script that uses
   `Bio.Entrez` has a line near the top:
   ```python
   Entrez.email = "your_email@example.com"  # required by NCBI; replace before running
   ```
   Replace the placeholder with your own email address before running.
2. Install ViennaRNA (RNAfold) and ensure it is on your system PATH. See
   https://www.tbi.univie.ac.at/RNA/ for installation instructions.
3. Download the following reference files into `feature_pipeline/data/`:
   - `SGD_all_ORFs_5prime_UTRs.fsa` and `SGD_all_ORFs_3prime_UTRs.fsa`
     (bulk UTR sequences, SGD downloads page)
   - `orf_trans_all.fasta` and `orf_coding_all.fasta` (bulk ORF sequences, SGD)
   - `NIHMS552811-supplement-1.xlsx` (mRNA isoform half-life data;
     Geisberg et al., 2014, *Cell*, Table S1)
4. The full target gene list (1710 systematic ORF names) is provided in
   `data/gene_list.txt`.

## Scripts and data sources

| # | Script | Features computed | Primary data source / tool |
|---|--------|--------------------|------------------------------|
| 01 | `01_cai_tai.R` | CAI, tAI | `seqinr` (CAI reference table), `tAI` package (dos Reis et al., 2004); tRNA gene copy numbers for *S. cerevisiae* from GtRNAdb (Chan & Lowe, 2016) |
| 02 | `02_sequence_extract.py` | Genomic/mRNA/protein length, mRNA MW, protein MW | SGD backend API |
| 03 | `03_utr_lengths.R` | 5' / 3' UTR length | SGD UTR FASTA files |
| 04 | `04_folding_energy.R` | UTR minimum free energy (MFE) | ViennaRNA (RNAfold) |
| 05 | `05_kozak_score.py` | Kozak strength, uORF count | NCBI Entrez; custom heuristic scoring function (not from a published model) |
| 06 | `06_cds_content.py` | CDS nucleotide/dinucleotide composition | NCBI Entrez |
| 07 | `07_utr_content.py` | UTR nucleotide/dinucleotide composition | SGD UTR FASTA files |
| 08 | `08_ensemble_free_energy.R` | Ensemble free energy, mean base-pair probability | ViennaRNA (partition function mode) |
| 09 | `09_are_score.py` | AU-rich element (ARE) score, ATTTA pentamer count | SGD 3' UTR FASTA file |
| 10 | `10_isoforms.py` | 5' / 3' UTR isoform counts | SGD UTR FASTA files |
| 11 | `11_motifs.py` | UTR sequence motif counts | SGD UTR FASTA files |
| 12 | `12_aug_uorf_context.py` | uORF / AUG sequence context | NCBI Entrez, reference genome (NC_0011xx accessions) |
| 13 | `13_stress_features.py` | Codon-usage-derived stress features | Kazusa Codon Usage Database (Nakamura et al., 2000) |
| 14 | `14_oxidation_data.py` | Cysteine content, oxidative potential | NCBI Entrez |
| 15 | `15_half_life.py` | mRNA / protein half-life | Geisberg et al. (2014) mRNA half-life dataset; SGD protein half-life field; heuristic fallback for genes not covered by either (flagged in output) |
| 16 | `16_protein_features_1.py` | Isoelectric point, instability index, aliphatic index | Biopython `ProteinAnalysis` (ProtParam) |
| 17 | `17_protein_features_2.py` | GRAVY, aromaticity, net charge | Biopython `ProteinAnalysis` (ProtParam) |
| 18 | `18_ramp_and_nend_features.py` | N-end rule stability score, translational ramp (5' AT-richness) | SGD bulk protein/CDS FASTA files; N-end rule (Bachmair, Finley & Varshavsky, 1986); translational ramp (Tuller et al., 2010) |

## Notes on data provenance

- Feature 05 (Kozak score) uses a custom heuristic, not a previously
  published quantitative model; treat as an approximate ranking rather than
  an absolute measure of initiation efficiency.
- Feature 15 (half-life) uses a three-tier hierarchy: published data first,
  SGD-curated data second, and a heuristic estimate only when neither is
  available. The output file's `mRNA_source` and `Protein_source` columns
  record which tier was used for each gene.
- All 18 scripts were run independently and their outputs merged into the
  final 115-column feature table (`feature_pipeline/output/all_features_1710genes.xlsx`).
