import os
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd

# Required by NCBI Entrez -- replace with your own email before running
Entrez.email = "your_email@example.com"

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def has_upstream_aug(utr_seq):
    return 'ATG' in utr_seq


def find_uorfs_with_kozak(utr_seq):
    """Find uAUGs with an in-frame stop codon, flagging those with a strong yeast Kozak context."""
    utr_seq = utr_seq.upper()
    stops = ['TAA', 'TAG', 'TGA']
    uorfs = []
    for i in range(len(utr_seq) - 2):
        if utr_seq[i:i + 3] == 'ATG':
            kozak_start = max(0, i - 10)
            kozak_seq = utr_seq[kozak_start:i]
            a_count = kozak_seq.count('A')
            g_count = kozak_seq.count('G')
            is_strong_kozak = (len(kozak_seq) >= 5) and (a_count >= 6) and (g_count <= 1)  # yeast-specific rule

            uorf_length = None
            for j in range(i + 3, len(utr_seq) - 2, 3):
                codon = utr_seq[j:j + 3]
                if len(codon) < 3:
                    break
                if codon in stops:
                    uorf_length = (j - i + 3) // 3
                    break

            if uorf_length:
                uorfs.append({
                    'position': i,
                    'strong_kozak': is_strong_kozak,
                    'length': uorf_length
                })
    return uorfs


with open(os.path.join(DATA_DIR, "gene_list.txt")) as f:
    orfs = [line.strip() for line in f if line.strip()]

results = []
for orf in orfs:
    try:
        term = f'{orf}[Gene Name] AND "Saccharomyces cerevisiae S288C"[Organism]'
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

        # 1000 bp upstream as a UTR proxy
        if strand_val == 'plus':
            utr_start = start - 1000
            utr_end = start - 1
            seq_strand = 1
        else:
            utr_start = end + 1
            utr_end = end + 1000
            seq_strand = 2  # reverse complement

        handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text",
                                strand=seq_strand, seq_start=utr_start, seq_stop=utr_end)
        fasta = handle.read()
        handle.close()
        utr_seq = ''.join(fasta.split('\n')[1:])

        aug_present = has_upstream_aug(utr_seq)
        uorfs = find_uorfs_with_kozak(utr_seq)
        uorf_present = len(uorfs) > 0
        strong_kozak_uorf = any(u['strong_kozak'] for u in uorfs)

        results.append({
            'ORF': orf,
            'Upstream AUG': aug_present,
            'uORF Present': uorf_present,
            'Strong Kozak uORF': strong_kozak_uorf
        })
        print(f"Processed {orf}")

    except Exception as e:
        print(f"Error for {orf}: {str(e)}")

df = pd.DataFrame(results)
df.to_excel(os.path.join(OUTPUT_DIR, "uorf_context.xlsx"), index=False)
print("Results saved to output/uorf_context.xlsx")
