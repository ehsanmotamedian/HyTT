import os
import time
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd

# Required by NCBI Entrez -- replace with your own email before running
Entrez.email = "your_email@example.com"

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

perc_lists = {p: [] for p in patterns}
orf_list = []

for orf in orfs:
    try:
        term = f'"{orf}"[Gene Name] AND "Saccharomyces cerevisiae S288C"[Organism]'
        handle = Entrez.esearch(db="gene", term=term)
        record = Entrez.read(handle)
        handle.close()

        if not record['IdList']:
            raise ValueError("No gene ID found")
        gene_id = record['IdList'][0]

        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        xml_data = handle.read()
        handle.close()
        tree = ET.fromstring(xml_data)

        accession = tree.find('.//Gene-commentary_accession').text
        from_pos = int(tree.find('.//Seq-interval_from').text)
        to_pos = int(tree.find('.//Seq-interval_to').text)
        strand_val = tree.find('.//Na-strand').get('value')

        start = from_pos + 1
        end = to_pos + 1

        if strand_val == 'plus':
            cds_start = start
            cds_end = end
            seq_strand = 1
        else:
            cds_start = end
            cds_end = start
            seq_strand = 2

        handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text",
                                strand=seq_strand, seq_start=cds_start, seq_stop=cds_end)
        fasta = handle.read()
        handle.close()
        cds_seq = ''.join(fasta.split('\n')[1:]).replace('\n', '')

        nuc_percs = calculate_nucleotide_percentages(cds_seq)
        for p in perc_lists:
            perc_lists[p].append(nuc_percs.get(p, 0.0))

        orf_list.append(orf)
        print(f"Processed {orf}")
        time.sleep(0.4)  # stay under NCBI's rate limit (3 requests/sec without an API key)

    except Exception as e:
        print(f"Error for {orf}: {str(e)}")
        orf_list.append(orf)
        for p in perc_lists:
            perc_lists[p].append(0.0)

data = {'ORF': orf_list}
for pattern in patterns:
    data[pattern] = perc_lists[pattern]
df = pd.DataFrame(data)
df.to_excel(os.path.join(OUTPUT_DIR, "cds_nucleotide_content.xlsx"), index=False)
print("Results saved to output/cds_nucleotide_content.xlsx")
print(f"Genes processed: {len(orf_list)}")
