import os
import io
import math
from collections import defaultdict
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
import xml.etree.ElementTree as ET
import pandas as pd

# Required by NCBI Entrez -- replace with your own email before running
Entrez.email = "your_email@example.com"

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Codon usage table for S. cerevisiae, Kazusa Codon Usage Database
# (Nakamura, Gojobori & Ikemura, 2000). Keys use RNA notation as published.
codon_usage_table = {
    'UUU': {'frequency_per_thousand': 26.1, 'number': 170666},
    'UCU': {'frequency_per_thousand': 23.5, 'number': 153557},
    'UAU': {'frequency_per_thousand': 18.8, 'number': 122728},
    'UGU': {'frequency_per_thousand': 8.1, 'number': 52903},
    'UUC': {'frequency_per_thousand': 18.4, 'number': 120510},
    'UCC': {'frequency_per_thousand': 14.2, 'number': 92923},
    'UAC': {'frequency_per_thousand': 14.8, 'number': 96596},
    'UGC': {'frequency_per_thousand': 4.8, 'number': 31095},
    'UUA': {'frequency_per_thousand': 26.2, 'number': 170884},
    'UCA': {'frequency_per_thousand': 18.7, 'number': 122028},
    'UAA': {'frequency_per_thousand': 1.1, 'number': 6913},
    'UGA': {'frequency_per_thousand': 0.7, 'number': 4447},
    'UUG': {'frequency_per_thousand': 27.2, 'number': 177573},
    'UCG': {'frequency_per_thousand': 8.6, 'number': 55951},
    'UAG': {'frequency_per_thousand': 0.5, 'number': 3312},
    'UGG': {'frequency_per_thousand': 10.4, 'number': 67789},
    'CUU': {'frequency_per_thousand': 12.3, 'number': 80076},
    'CCU': {'frequency_per_thousand': 13.5, 'number': 88263},
    'CAU': {'frequency_per_thousand': 13.6, 'number': 89007},
    'CGU': {'frequency_per_thousand': 6.4, 'number': 41791},
    'CUC': {'frequency_per_thousand': 5.4, 'number': 35545},
    'CCC': {'frequency_per_thousand': 6.8, 'number': 44309},
    'CAC': {'frequency_per_thousand': 7.8, 'number': 50785},
    'CGC': {'frequency_per_thousand': 2.6, 'number': 16993},
    'CUA': {'frequency_per_thousand': 13.4, 'number': 87619},
    'CCA': {'frequency_per_thousand': 18.3, 'number': 119641},
    'CAA': {'frequency_per_thousand': 27.3, 'number': 178251},
    'CGA': {'frequency_per_thousand': 3.0, 'number': 19562},
    'CUG': {'frequency_per_thousand': 10.5, 'number': 68494},
    'CCG': {'frequency_per_thousand': 5.3, 'number': 34597},
    'CAG': {'frequency_per_thousand': 12.1, 'number': 79121},
    'CGG': {'frequency_per_thousand': 1.7, 'number': 11351},
    'AUU': {'frequency_per_thousand': 30.1, 'number': 196893},
    'ACU': {'frequency_per_thousand': 20.3, 'number': 132522},
    'AAU': {'frequency_per_thousand': 35.7, 'number': 233124},
    'AGU': {'frequency_per_thousand': 14.2, 'number': 92466},
    'AUC': {'frequency_per_thousand': 17.2, 'number': 112176},
    'ACC': {'frequency_per_thousand': 12.7, 'number': 83207},
    'AAC': {'frequency_per_thousand': 24.8, 'number': 162199},
    'AGC': {'frequency_per_thousand': 9.8, 'number': 63726},
    'AUA': {'frequency_per_thousand': 17.8, 'number': 116254},
    'ACA': {'frequency_per_thousand': 17.8, 'number': 116084},
    'AAA': {'frequency_per_thousand': 41.9, 'number': 273618},
    'AGA': {'frequency_per_thousand': 21.3, 'number': 139081},
    'AUG': {'frequency_per_thousand': 20.9, 'number': 136805},
    'ACG': {'frequency_per_thousand': 8.0, 'number': 52045},
    'AAG': {'frequency_per_thousand': 30.8, 'number': 201361},
    'AGG': {'frequency_per_thousand': 9.2, 'number': 60289},
    'GUU': {'frequency_per_thousand': 22.1, 'number': 144243},
    'GCU': {'frequency_per_thousand': 21.2, 'number': 138358},
    'GAU': {'frequency_per_thousand': 37.6, 'number': 245641},
    'GGU': {'frequency_per_thousand': 23.9, 'number': 156109},
    'GUC': {'frequency_per_thousand': 11.8, 'number': 76947},
    'GCC': {'frequency_per_thousand': 12.6, 'number': 82357},
    'GAC': {'frequency_per_thousand': 20.2, 'number': 132048},
    'GGC': {'frequency_per_thousand': 9.8, 'number': 63903},
    'GUA': {'frequency_per_thousand': 11.8, 'number': 76927},
    'GCA': {'frequency_per_thousand': 16.2, 'number': 105910},
    'GAA': {'frequency_per_thousand': 45.6, 'number': 297944},
    'GGA': {'frequency_per_thousand': 10.9, 'number': 71216},
    'GUG': {'frequency_per_thousand': 10.8, 'number': 70337},
    'GCG': {'frequency_per_thousand': 6.2, 'number': 40358},
    'GAG': {'frequency_per_thousand': 19.2, 'number': 125717},
    'GGG': {'frequency_per_thousand': 6.0, 'number': 39359}
}

# Standard codon -> amino acid (DNA notation, matches CDS sequences fetched from GenBank)
codon_to_aa = {
    'AAA': 'K', 'AAC': 'N', 'AAG': 'K', 'AAT': 'N',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AGA': 'R', 'AGC': 'S', 'AGG': 'R', 'AGT': 'S',
    'ATA': 'I', 'ATC': 'I', 'ATG': 'M', 'ATT': 'I',
    'CAA': 'Q', 'CAC': 'H', 'CAG': 'Q', 'CAT': 'H',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'GAA': 'E', 'GAC': 'D', 'GAG': 'E', 'GAT': 'D',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'TAA': '*', 'TAC': 'Y', 'TAG': '*', 'TAT': 'Y',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TGA': '*', 'TGC': 'C', 'TGG': 'W', 'TGT': 'C',
    'TTA': 'L', 'TTC': 'F', 'TTG': 'L', 'TTT': 'F',
}

# Compute relative synonymous codon usage weights (w) for CAI.
# NOTE: codon_usage_table uses RNA notation (U) as published by Kazusa,
# while codon_to_aa (and the CDS sequences from GenBank) use DNA notation
# (T). Both are converted to DNA notation here before matching -- without
# this step, no codon containing U/T matches and only the 27 codons made
# solely of A/C/G are ever assigned a weight, silently corrupting the CAI
# calculation for every gene.
max_freq = defaultdict(float)
for codon, data in codon_usage_table.items():
    dna_codon = codon.replace('U', 'T')
    aa = codon_to_aa.get(dna_codon)
    if aa and aa != '*':
        freq = data['frequency_per_thousand']
        if freq > max_freq[aa]:
            max_freq[aa] = freq

w = {}
for codon, data in codon_usage_table.items():
    dna_codon = codon.replace('U', 'T')
    aa = codon_to_aa.get(dna_codon)
    if aa and aa != '*':
        freq = data['frequency_per_thousand']
        w[dna_codon] = freq / max_freq[aa] if max_freq[aa] > 0 else 0


def calculate_cai(cds_seq, w):
    cds_seq = cds_seq.upper()
    if len(cds_seq) % 3 != 0:
        return None
    log_w_sum = 0
    count = 0
    for i in range(0, len(cds_seq), 3):
        codon = cds_seq[i:i + 3]
        if codon in w and w[codon] > 0:
            log_w_sum += math.log(w[codon])
            count += 1
    if count == 0:
        return 0
    return math.exp(log_w_sum / count)


def proxy_delta_g(utr_seq):
    """DeltaG proxy using Tm_GC as a stability indicator (higher Tm ~ more stable, lower DeltaG)."""
    if len(utr_seq) == 0:
        return 0
    rna_seq = Seq(utr_seq).transcribe()
    tm = mt.Tm_GC(rna_seq, strict=False)
    return -tm


def count_cys(protein_seq):
    return protein_seq.count('C')


def disordered_proxy(protein_seq):
    """Fraction of disorder-promoting residues (P, G, S, E, K, Q, A)."""
    disordered_aa = set('PGSEKQA')
    count = sum(1 for aa in protein_seq if aa in disordered_aa)
    return count / len(protein_seq) if len(protein_seq) > 0 else 0


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
        seq_strand = 1 if strand_val == 'plus' else 2

        handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text",
                                strand=seq_strand, seq_start=start, seq_stop=end)
        gb_record = handle.read()
        handle.close()

        gb_parsed = next(SeqIO.parse(io.StringIO(gb_record), 'genbank'))

        cds_seq = ''
        protein_seq = ''
        for feature in gb_parsed.features:
            if feature.type == 'CDS':
                cds_seq = feature.extract(gb_parsed.seq)
                if 'translation' in feature.qualifiers:
                    protein_seq = feature.qualifiers['translation'][0]
                break

        # 1000 bp upstream as a UTR proxy
        if strand_val == 'plus':
            utr_start = start - 1000
            utr_end = start - 1
            utr_strand = 1
        else:
            utr_start = end + 1
            utr_end = end + 1000
            utr_strand = 2

        handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text",
                                strand=utr_strand, seq_start=utr_start, seq_stop=utr_end)
        fasta = handle.read()
        handle.close()
        utr_seq = ''.join(fasta.split('\n')[1:])

        cai = calculate_cai(str(cds_seq), w)
        delta_g_proxy = proxy_delta_g(utr_seq)
        cys_count = count_cys(protein_seq)
        disordered_frac = disordered_proxy(protein_seq)

        results.append({
            'ORF': orf,
            'CAI': cai,
            'DeltaG_proxy ( -Tm_GC )': delta_g_proxy,
            'Cys_count': cys_count,
            'Disordered_fraction': disordered_frac
        })
        print(f"Processed {orf}")
    except Exception as e:
        print(f"Error for {orf}: {str(e)}")

df = pd.DataFrame(results)
df.to_excel(os.path.join(OUTPUT_DIR, "stress_features.xlsx"), index=False)
print("Results saved to output/stress_features.xlsx")
