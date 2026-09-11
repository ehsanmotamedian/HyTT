import os
import re
from Bio import SeqIO
import pandas as pd

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

fasta_file = os.path.join(DATA_DIR, "SGD_all_ORFs_3prime_UTRs.fsa")
output_excel = os.path.join(OUTPUT_DIR, "are_scores.xlsx")

with open(os.path.join(DATA_DIR, "gene_list.txt")) as f:
    target_genes_ordered = [line.strip() for line in f if line.strip()]
target_genes = set(target_genes_ordered)


def count_overlapping_motifs(seq, pattern):
    return len(re.findall(f'(?={pattern})', seq))


results = []
found_genes = set()
missing_genes = target_genes.copy()

print("Scanning FASTA file for target genes...")
for record in SeqIO.parse(fasta_file, "fasta"):
    seq = str(record.seq).upper()
    length = len(seq)
    if length < 20:
        continue

    gene_match = re.search(r'(Y[A-Z][A-Z]\d{3}[A-Z](?:-[A-Z])?)', record.description)
    if not gene_match:
        continue
    gene_id = gene_match.group(1)

    if gene_id not in target_genes:
        continue
    found_genes.add(gene_id)
    missing_genes.discard(gene_id)

    # AU-rich element motifs (classic mRNA-decay signals; Chen & Shyu, 1995)
    pentamers = count_overlapping_motifs(seq, "ATTTA")
    clustered = count_overlapping_motifs(seq, "ATTTATTTA")
    t_stretches = len(re.findall(r'T{5,}', seq))
    nonamers = count_overlapping_motifs(seq, "[AT][AT]ATTTA[AT][AT]")
    at_density = (seq.count('A') + seq.count('T')) / length * 100 if length > 0 else 0
    are_score = pentamers * 5 + clustered * 15 + nonamers * 10 + t_stretches * 3
    norm_score = are_score * 1000 / length if length > 0 else 0

    results.append({
        "gene_id": gene_id,
        "length": length,
        "ATTTA_pentamers": pentamers,
        "clustered_classII": clustered,
        "nonamer_like": nonamers,
        "T_stretches_ge5": t_stretches,
        "AT_percent": round(at_density, 2),
        "ARE_score_raw": are_score,
        "ARE_score_per_kb": round(norm_score, 2)
    })

if not results:
    print("No target genes found -- check the FASTA file.")
else:
    df = pd.DataFrame(results)
    df_mean = df.groupby("gene_id").mean(numeric_only=True).reset_index()
    df_mean["gene_id"] = pd.Categorical(df_mean["gene_id"], categories=target_genes_ordered, ordered=True)
    df_mean = df_mean.sort_values("gene_id")
    df_mean.to_excel(output_excel, index=False)

    print(f"\nDone. Results for {len(found_genes)} genes written to:\n{output_excel}")
    print(f"Genes not found ({len(missing_genes)}): {', '.join(sorted(missing_genes))}")
