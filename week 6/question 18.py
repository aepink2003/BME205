import sys

class PeptideKmers:
    ''' Generate all possible k-mers of a given length from a peptide sequence '''
    def __init__(self, peptide, k):
        self.peptide = peptide
        self.k = k

    def generate_kmers(self):
        return [self.peptide[i:i+self.k] for i in range(len(self.peptide) - self.k + 1)]

class PeptideSpectrum:
    ''' Calculate the mass spectrum of a peptide '''
    mass_of_amino_acids = {
        'A': 71, 'G': 57, 'M': 131, 'S': 87, 'C': 103,
        'H': 137, 'N': 114, 'T': 101, 'D': 115, 'I': 113,
        'P': 97, 'V': 99, 'E': 129, 'K': 128, 'Q': 128,
        'W': 186, 'F': 147, 'L': 113, 'R': 156, 'Y': 163
    }

    def __init__(self, peptide):
        self.peptide = peptide

    def calculate_spectrum(self):
        subsets = [''] + [sub for k in range(1, len(self.peptide)+1) for sub in PeptideKmers(self.peptide, k).generate_kmers()]
        return sorted([sum(self.mass_of_amino_acids.get(aa, 0) for aa in subset) for subset in subsets])

class CyclopeptideSolver:
    ''' Solve the Cyclopeptide Sequencing Problem given a spectrum '''
    def __init__(self, spectrum):
        self.spectrum = spectrum

    def find_sequences(self):
        active_peptides = [[]]
        final_sequences = []

        while active_peptides:
            extended_peptides = [pep + [mass] for pep in active_peptides for mass in set(PeptideSpectrum.mass_of_amino_acids.values())]
            active_peptides = []

            for peptide in extended_peptides:
                if sum(peptide) == max(self.spectrum):
                    if PeptideSpectrum(peptide).calculate_spectrum() == self.spectrum:
                        final_sequences.append(peptide)
                elif all(mass in self.spectrum for mass in PeptideSpectrum(peptide).calculate_spectrum()):
                    active_peptides.append(peptide)

        return final_sequences

def main():
    input_spectrum = [int(mass) for mass in sys.stdin.readline().strip().split()]
    solver = CyclopeptideSolver(input_spectrum)
    for sequence in solver.find_sequences():
        print('-'.join(map(str, sequence)), end=' ')

if __name__ == "__main__":
    main()
