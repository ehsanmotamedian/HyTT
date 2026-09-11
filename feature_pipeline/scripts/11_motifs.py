import os
import re
from Bio import SeqIO
import pandas as pd

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

five_prime_file = os.path.join(DATA_DIR, "SGD_all_ORFs_5prime_UTRs.fsa")
three_prime_file = os.path.join(DATA_DIR, "SGD_all_ORFs_3prime_UTRs.fsa")
output_excel = os.path.join(OUTPUT_DIR, "motif_counts.xlsx")

with open(os.path.join(DATA_DIR, "gene_list.txt")) as f:
    target_genes_ordered = [line.strip() for line in f if line.strip()]
target_genes = set(target_genes_ordered)

motifs = {
    '5UTR': 'AUG',        # uORF start codon (RNA)
    '3UTR_ARE': 'AUUUA',   # AU-rich element
    '3UTR_POLYA': 'UAUAUA'  # polyadenylation signal
}


def count_motifs_in_fasta(file_path, motif, target_genes):
    results = {gene: 0 for gene in target_genes}
    found_genes = set()
    print(f"Scanning {os.path.basename(file_path)} for motif '{motif}'...")
    for record in SeqIO.parse(file_path, "fasta"):
        seq = str(record.seq).upper().replace('T', 'U')  # DNA -> RNA
        if len(seq) < 20:
            continue
        gene_match = re.search(r'(Y[A-Z][A-Z]\d{3}[A-Z](?:-[A-Z])?)', record.description)
        if not gene_match:
            continue
        gene_id = gene_match.group(1)
        if gene_id in target_genes:
            motif_count = len(re.findall(motif, seq))
            results[gene_id] += motif_count  # summed across isoforms
            found_genes.add(gene_id)
    missing = target_genes - found_genes
    print(f"Found: {len(found_genes)}, missing: {len(missing)}")
    return results, missing


five_utr_results, missing_5utr = count_motifs_in_fasta(five_prime_file, motifs['5UTR'], target_genes)
three_utr_are_results, missing_3utr_are = count_motifs_in_fasta(three_prime_file, motifs['3UTR_ARE'], target_genes)
three_utr_polya_results, missing_3utr_polya = count_motifs_in_fasta(three_prime_file, motifs['3UTR_POLYA'], target_genes)

results = []
for gene in target_genes_ordered:
    results.append({
        "gene_id": gene,
        "5UTR_AUG_count": five_utr_results.get(gene, 0),
        "3UTR_AUUUA_count": three_utr_are_results.get(gene, 0),
        "3UTR_UAUAUA_count": three_utr_polya_results.get(gene, 0)
    })
df = pd.DataFrame(results)

if df.empty:
    print("No data found -- check the FASTA files.")
else:
    df.to_excel(output_excel, index=False)
    print(f"\nDone. Results written to:\n{output_excel}")

    all_missing = missing_5utr | missing_3utr_are | missing_3utr_polya
    if all_missing:
        print(f"Genes missing from at least one file ({len(all_missing)}): {', '.join(sorted(all_missing))}")
