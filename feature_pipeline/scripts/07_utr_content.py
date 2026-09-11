import os
from Bio import SeqIO
import pandas as pd

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

patterns = ['GC', 'TT', 'GT', 'AT', 'GG', 'C', 'A', 'T', 'G', 'TC', 'AA', 'TA', 'CT', 'AC']


def calculate_nucleotide_percentages(seq):
    """Single- and di-nucleotide percentages (overlapping counts for dinucleotides)."""
    seq = seq.upper()
    length = len(seq)
    if length == 0:
        return {p: 0.0 for p in patterns}

    percentages = {}
    for p in patterns:
        if len(p) == 1:
            count = seq.count(p)
            percentages[p] = (count / length) * 100 if length > 0 else 0.0
        elif len(p) == 2:
            count = sum(1 for i in range(length - 1) if seq[i:i + 2] == p)
            percentages[p] = (count / (length - 1)) * 100 if length > 1 else 0.0
    return percentages


with open(os.path.join(DATA_DIR, "gene_list.txt")) as f:
    orfs = [line.strip() for line in f if line.strip()]

path_5utr = os.path.join(DATA_DIR, "SGD_all_ORFs_5prime_UTRs.fsa")
path_3utr = os.path.join(DATA_DIR, "SGD_all_ORFs_3prime_UTRs.fsa")

utr5_dict = {}
for record in SeqIO.parse(path_5utr, "fasta"):
    parts = record.id.split('_')
    if len(parts) > 4:
        orf = parts[4]
        utr5_dict[orf] = str(record.seq).replace('\n', '')

utr3_dict = {}
for record in SeqIO.parse(path_3utr, "fasta"):
    parts = record.id.split('_')
    if len(parts) > 4:
        orf = parts[4]
        utr3_dict[orf] = str(record.seq).replace('\n', '')

output_path = os.path.join(OUTPUT_DIR, "utr_nucleotide_content.xlsx")
with pd.ExcelWriter(output_path) as writer:
    for utr_type, utr_dict in [("5_UTR", utr5_dict), ("3_UTR", utr3_dict)]:
        data = {'ORF': []}
        for p in patterns:
            data[p] = []

        for orf in orfs:
            seq = utr_dict.get(orf, '')
            nuc_percs = calculate_nucleotide_percentages(seq)
            data['ORF'].append(orf)
            for p in patterns:
                data[p].append(nuc_percs.get(p, 0.0))

        df = pd.DataFrame(data)
        df.to_excel(writer, sheet_name=utr_type, index=False)
        print(f"Processed {utr_type} for {len(orfs)} genes")

print(f"Results saved to {output_path}")
print(f"Total genes processed: {len(orfs)}")
