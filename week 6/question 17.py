import sys

# Mapping of amino acids to their respective masses
amino_to_mass = {
    "G": 57, "A": 71, "S": 87, "P": 97, "V": 99,
    "T": 101, "C": 103, "I": 113, "L": 113, "N": 114,
    "D": 115, "K": 128, "Q": 128, "E": 129, "M": 131,
    "H": 137, "F": 147, "R": 156, "Y": 163, "W": 186
}

def peptide_mass_calc(peptide):
    ''' Calculate the total mass of a peptide '''
    return sum(amino_to_mass.get(amino, 0) for amino in peptide)

def linear_fragment(peptide):
    ''' Generates all linear fragments of the peptide '''
    fragments = []
    for i in range(1, len(peptide)):
        for j in range(len(peptide) - i + 1):
            fragments.append(peptide[j:j + i])
    return fragments

def cyclospectrum(peptide):
    ''' Calculate the cyclospectrum of the peptide '''
    spectrum = [0]
    linear_fragments = linear_fragment(peptide)
    
    # Add masses of linear fragments
    for fragment in linear_fragments:
        spectrum.append(peptide_mass_calc(fragment))
    
    # Add mass of the whole peptide
    spectrum.append(peptide_mass_calc(peptide))
    
    return sorted(spectrum)

def main():
    # Read the sequence from standard input
    sequence = sys.stdin.readline().strip()

    # Calculate and print the cyclospectrum
    spectrum = cyclospectrum(sequence)
    print(' '.join(map(str, spectrum)))

if __name__ == "__main__":
    main()
