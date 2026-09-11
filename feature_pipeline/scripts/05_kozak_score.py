import os
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd

# Required by NCBI Entrez -- replace with your own email before running
Entrez.email = "your_email@example.com"

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def kozak_strength(upstream_seq):
    """Yeast-specific Kozak strength (P_init) from the -10 to -1 upstream region.
    Heuristic scoring (not from a published quantitative model); capped at 0.9."""
    upstream_seq = upstream_seq.upper()
    if len(upstream_seq) < 5:
        return 0.1  # Minimum score for short UTRs to avoid 0
    a_count = upstream_seq.count('A')
    g_count = upstream_seq.count('G')
    base_score = a_count / len(upstream_seq)  # A-rich preference
    if g_count > 1:
        base_score *= 0.7  # Milder penalty for G
    if len(upstream_seq) >= 3 and upstream_seq[-3] == 'A':
        base_score += 0.1
    return min(0.9, max(0.1, base_score))


def find_uorfs_with_kozak(utr_seq):
    """Find uAUGs with an in-frame stop codon (candidate uORFs) and score each with kozak_strength."""
    utr_seq = utr_seq.upper()
    stops = ['TAA', 'TAG', 'TGA']
    uorfs = []
    for i in range(len(utr_seq) - 2):
        if utr_seq[i:i + 3] == 'ATG':
            kozak_start = max(0, i - 10)
            upstream_kozak = utr_seq[kozak_start:i]
            p_init = kozak_strength(upstream_kozak)

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
                    'p_init': p_init,
                    'length': uorf_length
                })
    return uorfs


with open(os.path.join(DATA_DIR, "gene_list.txt")) as f:
    orfs = [line.strip() for line in f if line.strip()]

results = []
utr_proxy_length = 200  # average yeast 5' UTR is ~100-200 bp

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

        # UTR proxy + 6 bp into CDS
        if strand_val == 'plus':
            fetch_start = start - utr_proxy_length
            fetch_end = start + 5
            seq_strand = 1
        else:
            fetch_start = end - 5
            fetch_end = end + utr_proxy_length
            seq_strand = 2  # reverse complement

        handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text",
                                strand=seq_strand, seq_start=fetch_start, seq_stop=fetch_end)
        fasta = handle.read()
        handle.close()
        full_seq = ''.join(fasta.split('\n')[1:])

        utr_seq = full_seq[:utr_proxy_length]
        aug = full_seq[utr_proxy_length:utr_proxy_length + 3]

        if aug != 'ATG':
            raise ValueError("Expected AUG not found at start position")

        main_upstream_kozak = utr_seq[-10:]
        main_p_init = kozak_strength(main_upstream_kozak)

        uorfs = find_uorfs_with_kozak(utr_seq)
        uorfs_sorted = sorted(uorfs, key=lambda x: x['position'])

        # P_skip via leaky scanning model: each uORF traps ribosomes with probability p_init
        p_skip = 1.0
        for u in uorfs_sorted:
            p_skip *= (1 - u['p_init'])

        efficiency = p_skip * main_p_init

        results.append({
            'ORF': orf,
            'Main Kozak Score': main_p_init,
            'Efficiency Estimate': efficiency,
            'uORFs': len(uorfs)
        })
        print(f"Processed {orf}")
    except Exception as e:
        print(f"Error for {orf}: {str(e)}")

df = pd.DataFrame(results)
df.to_excel(os.path.join(OUTPUT_DIR, "kozak_scores.xlsx"), index=False)
print("Results saved to output/kozak_scores.xlsx")
