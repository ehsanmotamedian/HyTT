import os
import requests
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def get_data(orf_name):
    url = f"https://www.yeastgenome.org/backend/locus/{orf_name}/sequence_details"
    response = requests.get(url)
    if response.status_code != 200:
        raise ValueError(f"Error fetching data for {orf_name}: {response.status_code}")
    data = response.json()

    genomic_length = None
    cds_residues = None
    protein_residues = None

    for seq in data.get('genomic_dna', []):
        if seq['strain']['display_name'] == 'S288C':
            genomic_length = len(seq['residues'])
            break
    for seq in data.get('coding_dna', []):
        if seq['strain']['display_name'] == 'S288C':
            cds_residues = seq['residues']
            break
    for seq in data.get('protein', []):
        if seq['strain']['display_name'] == 'S288C':
            protein_residues = seq['residues']
            break

    if genomic_length is None or cds_residues is None or protein_residues is None:
        raise ValueError(f"Data not found for S288C strain in {orf_name}")

    mrna_length = len(cds_residues)
    protein_length = len(protein_residues)

    # Strip trailing '*' from protein sequence if present
    protein_residues = protein_residues.rstrip('*')

    rna_seq = cds_residues.replace('T', 'U')
    mrna_mw = molecular_weight(Seq(rna_seq), seq_type="RNA")
    protein_mw = molecular_weight(Seq(protein_residues), seq_type="protein")

    return genomic_length, mrna_length, protein_length, mrna_mw, protein_mw


with open(os.path.join(DATA_DIR, "gene_list.txt")) as f:
    orf_names = [line.strip() for line in f if line.strip()]

data = []
for orf in orf_names:
    try:
        gene_len, mrna_len, prot_len, mrna_mw, prot_mw = get_data(orf)
        data.append({
            'ORF': orf,
            'Gene Length': gene_len,
            'mRNA Length': mrna_len,
            'Protein Length': prot_len,
            'mRNA MW (Da)': mrna_mw,
            'Protein MW (Da)': prot_mw
        })
        print(f"Processed: {orf}")
    except Exception as e:
        print(f"Error for {orf}: {str(e)}")

if data:
    df = pd.DataFrame(data)
    df.to_excel(os.path.join(OUTPUT_DIR, "sequence_data.xlsx"), index=False)
    print("Data saved to output/sequence_data.xlsx")
else:
    print("No data to save.")
