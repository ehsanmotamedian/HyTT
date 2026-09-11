from Bio import Entrez, SeqIO
from Bio.SeqUtils import ProtParam
from Bio.Data import CodonTable
import pandas as pd
import math
import json
import os
import time
from urllib.request import urlopen

# Required by NCBI Entrez -- replace with your own email before running
Entrez.email = "your_email@example.com"

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

with open(os.path.join(DATA_DIR, "gene_list.txt")) as f:
    genes = [line.strip() for line in f if line.strip()]

# Yeast codon usage table, Kazusa Codon Usage Database (frequency per 1000)
codon_table = {
    'TTT': 26.1, 'TTC': 18.4, 'TTA': 26.2, 'TTG': 27.2,
    'TCT': 23.5, 'TCC': 14.2, 'TCA': 18.7, 'TCG': 8.6,
    'TAT': 18.8, 'TAC': 14.8, 'TAA': 1.1, 'TAG': 0.5,
    'TGT': 8.1, 'TGC': 4.8, 'TGA': 0.7, 'TGG': 10.4,
    'CTT': 12.3, 'CTC': 5.4, 'CTA': 13.4, 'CTG': 10.5,
    'CCT': 13.5, 'CCC': 6.8, 'CCA': 18.3, 'CCG': 5.3,
    'CAT': 13.6, 'CAC': 7.8, 'CAA': 27.3, 'CAG': 12.1,
    'CGT': 6.4, 'CGC': 2.6, 'CGA': 3.0, 'CGG': 1.7,
    'ATT': 30.1, 'ATC': 17.2, 'ATA': 17.8, 'ATG': 20.9,
    'ACT': 20.3, 'ACC': 12.7, 'ACA': 17.8, 'ACG': 8.0,
    'AAT': 35.7, 'AAC': 24.8, 'AAA': 41.9, 'AAG': 30.8,
    'AGT': 14.2, 'AGC': 9.8, 'AGA': 21.3, 'AGG': 9.2,
    'GTT': 22.1, 'GTC': 11.8, 'GTA': 11.8, 'GTG': 10.8,
    'GCT': 21.2, 'GCC': 12.6, 'GCA': 16.2, 'GCG': 6.2,
    'GAT': 37.6, 'GAC': 20.2, 'GAA': 45.6, 'GAG': 19.2,
    'GGT': 23.9, 'GGC': 9.8, 'GGA': 10.9, 'GGG': 6.0
}

standard_table = CodonTable.unambiguous_dna_by_name["Standard"]

# Relative synonymous codon usage weights (w) from the frequency table
aa_to_codons = {}
for codon, aa in standard_table.forward_table.items():
    aa_to_codons.setdefault(aa, []).append(codon)
aa_to_codons['M'] = ['ATG']
aa_to_codons['W'] = ['TGG']

w = {}
for aa, codons_list in aa_to_codons.items():
    freqs = [codon_table.get(codon, 0) for codon in codons_list if codon in codon_table]
    max_f = max(freqs) if freqs and max(freqs) > 0 else 1.0
    for i, codon in enumerate(codons_list):
        if codon in codon_table:
            w[codon] = freqs[i] / max_f if max_f > 0 else 1.0


def calculate_cai(cds_seq, w):
    cds_seq = cds_seq.upper()
    if len(cds_seq) % 3 != 0 or len(cds_seq) < 3:
        return 0.0
    codons = [cds_seq[i:i + 3] for i in range(0, len(cds_seq), 3)]
    valid_codons = [c for c in codons if c in w and w[c] > 0]
    if not valid_codons:
        return 0.0
    sum_log = sum(math.log(w[c]) for c in valid_codons if w[c] > 0)
    return math.exp(sum_log / len(valid_codons))


# Primary mRNA half-life source: Geisberg et al. (2014), Cell, Table S1
# (download the supplementary file and place it in DATA_DIR).
mrna_supp_file = os.path.join(DATA_DIR, "NIHMS552811-supplement-1.xlsx")
try:
    mrna_data = pd.read_excel(mrna_supp_file)
    mrna_dict = dict(zip(mrna_data['Gene'], mrna_data['Half-life (min)']))
except FileNotFoundError:
    mrna_dict = {}
    print(f"mRNA half-life file not found at {mrna_supp_file} -- falling back to estimates for all genes.")

mRNA_half_lives = []
protein_half_lives = []
sources_mRNA = []
sources_protein = []

for gene in genes:
    # Tier 1: published mRNA half-life (Geisberg et al., 2014)
    mRNA_hl = mrna_dict.get(gene, None)
    source_mRNA = "Literature" if mRNA_hl is not None else "Estimated"
    prot_hl = None
    source_protein = "Estimated"

    # Tier 1: SGD-curated protein half-life
    url = f"https://www.yeastgenome.org/backend/locus/{gene}"
    try:
        with urlopen(url) as response:
            data = json.loads(response.read().decode())

        if 'protein_overview' in data and 'half_life' in data['protein_overview']:
            half_life_dict = data['protein_overview']['half_life']
            prot_hl_str = half_life_dict.get('data_value', '')
            if prot_hl_str:
                prot_hl_str_clean = prot_hl_str.replace('>=', '').replace(' ', '')
                try:
                    prot_hl = float(prot_hl_str_clean)
                    source_protein = "SGD"
                except ValueError:
                    print(f"Invalid SGD protein half-life value for {gene}: {prot_hl_str}")
    except Exception as e:
        print(f"Error fetching SGD data for {gene}: {e}")

    time.sleep(0.5)  # stay under SGD's rate limit

    # Tier 2 (fallback): heuristic estimates for genes not covered by the
    # primary sources above. NOTE: these are approximations devised for
    # pipeline completeness, not derived from a validated published model.
    try:
        handle = Entrez.esearch(db="gene", term=f"{gene}[Gene Name] AND Saccharomyces cerevisiae[Organism]")
        record = Entrez.read(handle)
        handle.close()
        time.sleep(0.3)
        if record["IdList"]:
            gene_id = record["IdList"][0]

            if mRNA_hl is None:
                handle = Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id, linkname="gene_nuccore_refseqrna")
                link_record = Entrez.read(handle)
                handle.close()
                time.sleep(0.3)
                if link_record and link_record[0].get("LinkSetDb"):
                    links = link_record[0]["LinkSetDb"][0]["Link"]
                    if links:
                        mrna_id = links[0]["Id"]
                        handle = Entrez.efetch(db="nuccore", id=mrna_id, rettype="fasta_cds_na", retmode="text")
                        cds_seq = str(SeqIO.read(handle, "fasta").seq)
                        handle.close()
                        cai_value = calculate_cai(cds_seq, w)
                        # Linear interpolation of mRNA half-life from CAI
                        if cai_value < 0.4:
                            mRNA_hl = 5.0
                        elif cai_value > 0.7:
                            mRNA_hl = 18.0
                        else:
                            mRNA_hl = 5 + 13 * ((cai_value - 0.4) / 0.3)

            if prot_hl is None:
                handle = Entrez.elink(dbfrom="gene", db="protein", id=gene_id, linkname="gene_protein_refseq")
                link_record = Entrez.read(handle)
                handle.close()
                time.sleep(0.3)
                if link_record and link_record[0].get("LinkSetDb"):
                    links = link_record[0]["LinkSetDb"][0]["Link"]
                    if links:
                        prot_id = links[0]["Id"]
                        handle = Entrez.efetch(db="protein", id=prot_id, rettype="fasta", retmode="text")
                        prot_seq = str(SeqIO.read(handle, "fasta").seq)
                        handle.close()
                        analysis = ProtParam.ProteinAnalysis(prot_seq)
                        index = analysis.instability_index()
                        # Empirical decay relationship from instability index
                        prot_hl = 20 * math.exp(-0.05 * index)
    except Exception as e:
        print(f"Error estimating half-life for {gene}: {e}")

    time.sleep(0.5)

    mRNA_half_lives.append(mRNA_hl)
    protein_half_lives.append(prot_hl)
    sources_mRNA.append(source_mRNA)
    sources_protein.append(source_protein)

df = pd.DataFrame({
    "Gene": genes,
    "mRNA_half_life_min": mRNA_half_lives,
    "Protein_half_life_hours": protein_half_lives,
    "mRNA_source": sources_mRNA,
    "Protein_source": sources_protein
})
df.to_excel(os.path.join(OUTPUT_DIR, "gene_half_lives.xlsx"), index=False)
print("Results saved to output/gene_half_lives.xlsx")
