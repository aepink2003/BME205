import math
from random import randint
import sys


class FastAreader:
    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''

        self.fname = fname
        self.fileH = None

    def doOpen(self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):

        header = ''
        sequence = ''

        with self.doOpen() as self.fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = self.fileH.readline()
            while not line.startswith('>'):
                line = self.fileH.readline()
            header = line[1:].rstrip()

            for line in self.fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header, sequence



class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='commandLine BME205 solution for the Finding CRISPR arrays assignments',
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )

        self.parser.add_argument('-i', '--iterations', type=int, required=True, action='store', help='Number of iterations')
        self.parser.add_argument('-k', '--motifLength', type=int, required=True, action='store', help='Motif length')
        self.parser.add_argument('-p', '--pseudocount', type=float, required=True, action='store', help='Pseudocount value')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)



class Usage(Exception):
    '''
    Signal a Usage error, evoking a usage statement and eventual exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg



class Motif():
    '''
    '''
    def __init__(self, sequences, iterations, motif_length, pseudocount):
        '''
        '''
        self.sequences = sequences
        self.iterations = int(iterations)
        self.motif_length = int(motif_length)
        self.pseudocount = int(pseudocount)

    
    def generate_profiles_with_motif_and_pseudocounts(self, DNA_sequences, pseudocount, motif_strings):
        """
        Generates profiles with given motif and pseudocounts.

        Args:
            DNA_sequences (list): A list of DNA sequences.
            pseudocount (int): Integer value for pseudocounts.
            motif_strings (list): List of strings of initiated motifs.

        Returns:
            dict: Dictionary with 'A', 'C', 'G', 'T' as keys.
        """

        # Calculate the total number of sequences after adding pseudocounts
        n = 4 * pseudocount + len(DNA_sequences)

        # Initialize a dictionary to store the profile
        profile = {'A': [], 'C': [], 'G': [], 'T': []}

        # Iterate through positions in the motifs
        for i in range(len(motif_strings[0])):
            # Extract the i-th position from all motifs
            temp = [motif[i] for motif in motif_strings]

            # Calculate probabilities with pseudocounts
            profile['A'].append((temp.count('A') + pseudocount) / n)
            profile['C'].append((temp.count('C') + pseudocount) / n)
            profile['G'].append((temp.count('G') + pseudocount) / n)
            profile['T'].append((temp.count('T') + pseudocount) / n)

            # Clear the temporary list
            temp.clear()

        return profile

    

    def get_new_motif_based_on_profile(self, DNA_sequences, profile, k):
        """
        Gets a new motif based on a given profile.

        Args:
            DNA_sequences (list): A list of DNA sequences.
            profile (dict): Profile in dictionary form.
            k (int): Integer representing k-mer.

        Returns:
            list: List of strings representing motifs.
        """

        # Initialize an empty list to store motifs
        motifs = []

        # Iterate through each DNA sequence
        for seq in DNA_sequences:
            sub_sequences = []

            # Generate all possible k-mers in the sequence
            for x in range(0, len(seq)-k+1):
                sub_sequences.append(seq[x:x+k])

            # Calculate probability of each possible motif
            probabilities = []
            for i in range(len(sub_sequences)):
                sub_seq = sub_sequences[i]
                prob_list = [profile[sub_seq[j]][j] for j in range(k)]
                probabilities.append(math.prod(prob_list))

            # Find motif with maximum probability
            max_probability = max(probabilities)
            max_index = probabilities.index(max_probability)
            best_motif = sub_sequences[max_index]

            # Append the best motif to the list
            motifs.append(best_motif)

        return motifs

    

    
    def calculate_relative_entropy_score(self, profile, DNA_sequences, pseudocount, k):
        """
        Uses relative entropy method to get relative distance with null model.

        Args:
            profile (dict): Dictionary form of profile.
            DNA_sequences (list): A list of DNA sequences.
            pseudocount (int): Integer value for pseudocounts.
            k (int): Integer representing k-mer length.

        Returns:
            float: Score based on relative entropy.
        """

        # Construct null model
        concatenated_sequence = ''.join(DNA_sequences)
        null_model = {}  # Q(i)
        total_length = len(concatenated_sequence)
        null_model['A'] = (concatenated_sequence.count('A') + pseudocount) / (total_length + 4 * pseudocount)
        null_model['C'] = (concatenated_sequence.count('C') + pseudocount) / (total_length + 4 * pseudocount)
        null_model['G'] = (concatenated_sequence.count('G') + pseudocount) / (total_length + 4 * pseudocount)
        null_model['T'] = (concatenated_sequence.count('T') + pseudocount) / (total_length + 4 * pseudocount)

        # Calculate relative entropy score
        score = 0
        for i in range(k):
            for base in list(profile.keys()):
                score += profile[base][i] * math.log2(profile[base][i] / null_model[base])

        return score
            

    def randomized_motif_search(self, DNA_sequences, k, pseudocount):
        """
        Perform Randomized Motif Search algorithm.

        Args:
            DNA_sequences (list): A list of DNA sequences.
            k (int): Integer representing k-mer length.
            pseudocount (int): Integer value for pseudocounts.

        Returns:
            list: List of best motifs.
        """        

        # Initialize motifs with random selections
        motifs = []
        for sequence in DNA_sequences:
            start_position = randint(0, len(sequence) - k)
            motifs.append(sequence[start_position:start_position + k])

        best_score = 0
        best_motifs = motifs

        while True:
            # Calculate current profile
            current_profile = self.generate_profiles_with_motif_and_pseudocounts(DNA_sequences, pseudocount, motifs)

            # Get new motifs based on current profile
            motifs = self.get_new_motif_based_on_profile(DNA_sequences, current_profile, k)

            # Calculate current score
            current_score = self.calculate_relative_entropy_score(self.generate_profiles_with_motif_and_pseudocounts(DNA_sequences, pseudocount, motifs), DNA_sequences, pseudocount, k)

            # Update best score and motifs if necessary
            if current_score > best_score:
                best_score = current_score
                best_motifs = motifs
            else:
                return best_motifs




     
def main (inFile=None, options = None):

    """
    Main function for executing the motif finding algorithm.

    Args:
        inFile (str): Path to the input file.
        options (obj): Options object for command-line arguments.

    Returns:
        None
    """
    
    # initialize command line object
    cl = CommandLine(options)

    if sys.stdin:

        # Array to store all the sequences in the fasta file.
        sequences = []

        # Read all sequences and store then in the 'sequences' list
        for header, seq in FastAreader().readFasta():
            sequences.append(seq)
                
        # Initialize motif object
        motif = Motif(sequences, cl.args.iterations, cl.args.motifLength, cl.args.pseudocount)
        
        # Iterate
        best_motifs = []
        best_score = 0
        for _ in range(cl.args.iterations):
            current_motif = motif.randomized_motif_search(sequences, cl.args.motifLength, cl.args.pseudocount)
            current_profile = motif.generate_profiles_with_motif_and_pseudocounts(sequences, cl.args.pseudocount, current_motif)
            current_score = motif.calculate_relative_entropy_score(current_profile, sequences, cl.args.pseudocount, cl.args.motifLength)
            if current_score > best_score:
                best_score = current_score
                best_motifs = current_motif

        # Find consensus
        temp_profile = motif.generate_profiles_with_motif_and_pseudocounts(sequences, cl.args.pseudocount, best_motifs)
        consensus = ''
        for i in range(0, cl.args.motifLength):
            temp = []
            for k in list(temp_profile.keys()):
                temp.append(temp_profile[k][i])
            max_value = max(temp)
            max_index = temp.index(max_value)
            if max_index == 0:
                consensus += 'A'
            elif max_index == 1:
                consensus += 'C'
            elif max_index == 2:
                consensus += 'G'
            else:
                consensus += 'T'

        
        print(inFile + ': ')
        print(best_score)
        print(consensus + '\n')

    else:
        raise Usage("Please enter the path of the STDIN file: randomizedMotifSearch.py -i 1000 -k 13 -p 1 < input.fa > output.out")
    
       
if __name__ == "__main__":
    main(inFile='p1860Crisprs', options=['-i 1000', '-k 13', '-p 1'])
    

