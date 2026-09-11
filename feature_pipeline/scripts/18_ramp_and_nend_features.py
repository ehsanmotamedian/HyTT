import os
import gzip
import pandas as pd
from Bio import SeqIO

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# SGD bulk protein/CDS FASTA files (place in DATA_DIR; .gz is auto-detected)
protein_file = os.path.join(DATA_DIR, "orf_trans_all.fasta")
cds_file = os.path.join(DATA_DIR, "orf_coding_all.fasta")

with open(os.path.join(DATA_DIR, "gene_list.txt")) as f:
    genes = [line.strip() for line in f if line.strip()]


def open_file(path):
    if path.endswith('.gz'):
        return gzip.open(path, "rt")
    return open(path, "r")


prot_dict = {}
with open_file(protein_file) as f:
    for record in SeqIO.parse(f, "fasta"):
        gene = record.id.split('|')[0].lstrip('>')
        prot_dict[gene] = str(record.seq)

cds_dict = {}
with open_file(cds_file) as f:
    for record in SeqIO.parse(f, "fasta"):
        gene = record.id.split('|')[0].lstrip('>')
        cds_dict[gene] = str(record.seq).upper()

results = []
for gene in genes:
    seq_prot = prot_dict.get(gene)
    cds_seq = cds_dict.get(gene)

    first_aa = seq_prot[0] if seq_prot and len(seq_prot) > 0 else None

    # N-end rule stability score (Bachmair, Finley & Varshavsky, 1986):
    # stabilizing, destabilizing, and intermediate N-terminal residue classes
    if first_aa in 'MAVSTG':
        n_score = 1.0
    elif first_aa in 'RKHLFWI':
        n_score = -1.0
    elif first_aa in 'DENQC':
        n_score = -0.5
    else:
        n_score = 0.0 if first_aa else None

    # Translational ramp (Tuller et al., 2010): AT-richness of the first 50
    # codons relative to the rest of the CDS, associated with a slow-
    # translating "ramp" that reduces ribosome collisions downstream
    at_first50 = None
    ramp_score = None
    if cds_seq and len(cds_seq) >= 150:
        first_150 = cds_seq[:150]
        rest = cds_seq[150:]
        at_first = (first_150.count('A') + first_150.count('T')) / 150.0
        at_rest = (rest.count('A') + rest.count('T')) / len(rest) if rest else at_first
        at_first50 = round(at_first * 100, 2)
        ramp_score = round(at_first - at_rest, 4)

    results.append({
        "Gene": gene,
        "First_AA": first_aa,
        "N_end_stability_score": n_score,
        "AT_content_first_50_codons_%": at_first50,
        "Ramp_score": ramp_score,
        "Protein_length": len(seq_prot) if seq_prot else None,
        "CDS_length": len(cds_seq) if cds_seq else None
    })

df = pd.DataFrame(results)
df.to_excel(os.path.join(OUTPUT_DIR, "ramp_and_nend_features.xlsx"), index=False)
print("Results saved to output/ramp_and_nend_features.xlsx")
