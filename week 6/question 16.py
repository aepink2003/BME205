import sys

# Restructured DNA codons dictionary
codon_pairs = [ 
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATA": "I", "ATT": "I", "ATC": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "_", "TAG": "_", "TGA": "_"
]
DNA_Codons = dict(codon_pairs)
inverse_Codons = {value: key for key, value in DNA_Codons.items()}

def complement_dna(sequence):
    translation_table = str.maketrans('ACGTN', 'TGCAN')
    return sequence[::-1].translate(translation_table)

def split_into_codons(dna_sequence, offset=0):
    return [dna_sequence[i:i + 3] for i in range(offset, len(dna_sequence), 3)]

def translate_to_amino_acids(dna_seq):
    return ''.join(DNA_Codons.get(codon, '') for codon in split_into_codons(dna_seq))

def search_amino_acid_sequences(dna, amino_acid_seq):
    found_sequences = []
    dna = dna.upper()
    reversed_dna = complement_dna(dna)
    amino_acid_seq = amino_acid_seq.upper()
    length = len(amino_acid_seq) * 3

    for i in range(len(dna) - length + 1):
        segment = dna[i:i + length]
        if translate_to_amino_acids(segment) == amino_acid_seq:
            found_sequences.append(segment)
    
    for i in range(len(reversed_dna) - length + 1):
        segment = reversed_dna[i:i + length]
        if translate_to_amino_acids(segment) == amino_acid_seq:
            found_sequences.append(complement_dna(segment))

    return found_sequences

def run_amino_acid_search():
    dna_sequence, amino_acid_sequence = (line.strip().upper() for line in sys.stdin)
    matching_sequences = search_amino_acid_sequences(dna_sequence, amino_acid_sequence)

    for sequence in matching_sequences:
        print(sequence)

if __name__ == "__main__":
    run_amino_acid_search()
