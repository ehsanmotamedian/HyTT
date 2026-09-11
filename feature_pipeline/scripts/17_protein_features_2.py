import os
import time
import numpy as np
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Required by NCBI Entrez -- replace with your own email before running
Entrez.email = "your_email@example.com"

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

with open(os.path.join(DATA_DIR, "gene_list.txt")) as f:
    genes = [line.strip() for line in f if line.strip()]


def fetch_protein_sequence(gene):
    try:
        handle = Entrez.esearch(db="protein", term=f'"{gene}"[Gene] AND Saccharomyces cerevisiae[Organism]', retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if record['IdList']:
            prot_id = record['IdList'][0]
            handle = Entrez.efetch(db="protein", id=prot_id, rettype="fasta", retmode="text")
            seq = str(SeqIO.read(handle, "fasta").seq)
            handle.close()
            return seq
    except Exception:
        return None
    return None


def calculate_aliphatic_index(prot):
    """Aliphatic index (Ikai, 1980): relative volume occupied by aliphatic side chains."""
    percents = prot.amino_acids_percent
    x_a = percents.get('A', 0) * 100
    x_v = percents.get('V', 0) * 100
    x_i = percents.get('I', 0) * 100
    x_l = percents.get('L', 0) * 100
    return x_a + 2.9 * x_v + 3.9 * (x_i + x_l)


data = []
for gene in genes:
    time.sleep(0.34)  # stay under NCBI's rate limit
    prot_seq = fetch_protein_sequence(gene)

    if prot_seq:
        prot = ProteinAnalysis(prot_seq)
        row = {
            'Gene': gene,
            'Isoelectric_point_pI': round(prot.isoelectric_point(), 3),
            'Instability_index': round(prot.instability_index(), 3),
            'Aliphatic_index': round(calculate_aliphatic_index(prot), 3),
            'Protein_length': len(prot_seq)
        }
        data.append(row)
    else:
        print(f"Failed: {gene}")
        data.append({'Gene': gene, 'Isoelectric_point_pI': np.nan, 'Instability_index': np.nan,
                      'Aliphatic_index': np.nan, 'Protein_length': np.nan})

df = pd.DataFrame(data)
df.to_excel(os.path.join(OUTPUT_DIR, "protein_features_2.xlsx"), index=False)
print("Results saved to output/protein_features_2.xlsx")
print(f"Genes with a resolved protein sequence: {sum(df['Isoelectric_point_pI'].notna())} of {len(genes)}")
