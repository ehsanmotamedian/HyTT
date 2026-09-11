import os
import io
from Bio import Entrez, SeqIO
import xml.etree.ElementTree as ET
import pandas as pd

# Required by NCBI Entrez -- replace with your own email before running
Entrez.email = "your_email@example.com"

DATA_DIR = "data"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

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
                cds_seq = str(feature.extract(gb_parsed.seq))
                if 'translation' in feature.qualifiers:
                    protein_seq = feature.qualifiers['translation'][0]
                break

        # mRNA oxidative potential proxy: G content (guanine is the most
        # oxidation-susceptible base, e.g. via 8-oxo-guanine formation)
        if cds_seq:
            g_count_mrna = cds_seq.upper().count('G')
            cds_length = len(cds_seq)
            mrna_oxidative_potential = g_count_mrna / cds_length if cds_length > 0 else 0
        else:
            mrna_oxidative_potential = 0

        # Protein oxidative potential proxy: cysteine content (most
        # oxidation-susceptible residue)
        cys_count_protein = protein_seq.count('C')
        protein_length = len(protein_seq)
        protein_oxidative_potential = cys_count_protein / protein_length if protein_length > 0 else 0

        oxidative_ratio = protein_oxidative_potential / mrna_oxidative_potential if mrna_oxidative_potential > 0 else 0

        results.append({
            'ORF': orf,
            'mRNA_Oxidative_Potential': mrna_oxidative_potential,
            'Protein_Oxidative_Potential': protein_oxidative_potential,
            'Oxidative_Ratio': oxidative_ratio
        })
        print(f"Processed {orf}")

    except Exception as e:
        print(f"Error for {orf}: {str(e)}")

df = pd.DataFrame(results)
df.to_excel(os.path.join(OUTPUT_DIR, "oxidative_potential.xlsx"), index=False)
print("Results saved to output/oxidative_potential.xlsx")
