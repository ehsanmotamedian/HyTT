import os
import re
from Bio import SeqIO
import pandas as pd

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

fasta_3utr = os.path.join(DATA_DIR, "SGD_all_ORFs_3prime_UTRs.fsa")
fasta_5utr = os.path.join(DATA_DIR, "SGD_all_ORFs_5prime_UTRs.fsa")
output_excel = os.path.join(OUTPUT_DIR, "utr_isoform_counts.xlsx")

with open(os.path.join(DATA_DIR, "gene_list.txt")) as f:
    target_genes_ordered = [line.strip() for line in f if line.strip()]
target_genes = set(target_genes_ordered)


def extract_gene_id(description):
    gene_match = re.search(r'(Y[A-Z][A-Z]\d{3}[A-Z](?:-[A-Z])?)', description)
    return gene_match.group(1) if gene_match else None


def count_isoforms(fasta_file, target_genes):
    counts = {gene: 0 for gene in target_genes}
    found_genes = set()
    print(f"Scanning {os.path.basename(fasta_file)} for isoforms...")
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        if len(seq) < 20:
            continue
        gene_id = extract_gene_id(record.description)
        if gene_id and gene_id in target_genes:
            counts[gene_id] += 1
            found_genes.add(gene_id)
    missing = target_genes - found_genes
    print(f"Found: {len(found_genes)}, missing: {len(missing)}")
    return counts, missing


counts_3utr, missing_3utr = count_isoforms(fasta_3utr, target_genes)
counts_5utr, missing_5utr = count_isoforms(fasta_5utr, target_genes)

results = []
for gene in target_genes_ordered:
    results.append({
        "gene_id": gene,
        "num_5UTR_isoforms": counts_5utr.get(gene, 0),
        "num_3UTR_isoforms": counts_3utr.get(gene, 0)
    })
df = pd.DataFrame(results)

if df.empty:
    print("No data found -- check the FASTA files.")
else:
    df.to_excel(output_excel, index=False)
    print(f"\nDone. Results written to:\n{output_excel}")

    all_missing = missing_3utr | missing_5utr
    if all_missing:
        print(f"Genes missing from at least one file ({len(all_missing)}): {', '.join(sorted(all_missing))}")
