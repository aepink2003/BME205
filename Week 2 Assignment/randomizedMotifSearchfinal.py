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
        self.parser.add_argument('-p', '--psuedocount_value', type=float, required=True, action='store', help='psuedocount_value value')

        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)



class Usage(Exception):
    '''
        tells user when something has gone terribly wrong unable to exec
    '''
    def __init__(self, msg):
        self.msg = msg



class Motif():
    '''
    Class to perform motif search in DNA sequences using the Randomized Motif Search algorithm.
    '''
    def __init__(self, seq, iterations, kmer_length, psuedocount_value):
        '''
        '''
        self.sequences = seq
        self.iterations = int(iterations)
        self.kmer_length = int(kmer_length)
        self.psuedocount_value = int(psuedocount_value)


    
    #idkkk what im doing here
    def rando_motif_search(self, DNA_seq, k, psuedocount_value):
        """
        Random Motif Search algo
        """        
        # init  random selections
        motif = []
        for seq in DNA_seq:
            a_pos = randint(0, len(seq) - k)
            motif.append(seq[a_pos:a_pos + k])
        w_score = 0
        w_motif = motif
        while True:
            prof = self.make_motif_psuedocount_prof(DNA_seq, psuedocount_value, motif)
            # get motifs based on current profile
            motif = self.newMotif_prof(DNA_seq, prof, k)
            # Calculate current score
            score = self.calc_rel_entropy(self.make_motif_psuedocount_prof(DNA_seq, psuedocount_value, motif), DNA_seq, psuedocount_value, k)
            # Update best score and motifs if necessary
            if score > w_score:
                w_score = score
                w_motif = motif
            else:
                return w_motif
    

    def newMotif_prof(self, DNA_seq, profile, k):
        """
        creates/gives a new motif based on a given profile
        """
        # init an empty list to store motifs
        motif = []
        for seq in DNA_seq:
            lower_seq = []
            # generate all possible k-mers in the sequence
            for x in range(0, len(seq)-k+1):
                lower_seq.append(seq[x:x+k])
            # calc probability of possible motif
            prob = []
            for i in range(len(lower_seq)):
                lower_seq = lower_seq[i]
                prob_list = [profile[lower_seq[j]][j] for j in range(k)]
                prob.append(math.prod(prob_list))
            # maximum prob
            max_prob = max(prob)
            maxi = prob.index(max_prob)
            best_motif = lower_seq[maxi]
            # Append to list
            motif.append(best_motif)

        return motif

    

    
    def calc_rel_entropy(self, motif_profile, DNA_seq, psuedocount_value, k):
        """
        calcs rel. entropy usiing relative entropy  

        """
        # KK helped me with this
        concat_seq = ''.join(DNA_seq)
        null_model = {} 
        totalLen = len(concat_seq)
        null_model['A'] = (concat_seq.count('A') + psuedocount_value) / (totalLen + 4 * psuedocount_value)
        null_model['C'] = (concat_seq.count('C') + psuedocount_value) / (totalLen + 4 * psuedocount_value)
        null_model['G'] = (concat_seq.count('G') + psuedocount_value) / (totalLen + 4 * psuedocount_value)
        null_model['T'] = (concat_seq.count('T') + psuedocount_value) / (totalLen + 4 * psuedocount_value)
        #relative entropy score
        ent_score = 0
        for i in range(k):
            for base in list(motif_profile.keys()):
                score += motif_profile[base][i] * math.log2(motif_profile[base][i] / null_model[base])

        return ent_score
    
    def make_motif_psuedocount_prof(self, DNA_seq, psuedocount_value, motif_strings):
        """
        makes profile with motif and pc value
        """
        # calc the total              sequences despues psuedocount_values
        s = 4 * psuedocount_value + len(DNA_seq)
        # init dictionary 
        profile = {'A': [], 'C': [], 'G': [], 'T': []}
        # iterate through pos. in the motifs
        for i in range(len(motif_strings[0])):
            # take  i-th position from all motifs
            hehelist = [motif[i] for motif in motif_strings]
            # calc prob 
            profile['A'].append((hehelist.count('A') + psuedocount_value) / s)
            profile['C'].append((hehelist.count('C') + psuedocount_value) / s)
            profile['G'].append((hehelist.count('G') + psuedocount_value) / s)
            profile['T'].append((hehelist.count('T') + psuedocount_value) / s)

            # clear hehelist list
            hehelist.clear()

        return profile

            







     
def main (inFile=None, options = None):

    """
    Main function to execute this whole thing
    """
    cmdl = CommandLine(options)
    if sys.stdin:
        #store all the sequences in the fasta file.
        seq = []
        # store then in seq list
        for header, seq in FastAreader().readFasta():
            seq.append(seq)               
        # init motif object
        motif = Motif(seq, cmdl.args.iterations, cmdl.args.motifLength, cmdl.args.psuedocount_value)
        w_motif = []
        w_score = 0
        for _ in range(cmdl.args.iterations):
            current_motif = motif.rando_motif_search(seq, cmdl.args.motifLength, cmdl.args.psuedocount_value)
            prof = motif.make_motif_psuedocount_prof(seq, cmdl.args.psuedocount_value, current_motif)
            score = motif.calc_rel_entropy(prof, seq, cmdl.args.psuedocount_value, cmdl.args.motifLength)
            if score > w_score:
                w_score = score
                w_motif = current_motif
        #consensus
        tt_prof = motif.make_motif_psuedocount_prof(seq, cmdl.args.psuedocount_value, w_motif)
        consensus = ''
        for i in range(0, cmdl.args.motifLength):
            hehelist = []
            for k in list(tt_prof.keys()):
                hehelist.append(tt_prof[k][i])
            max_val = max(hehelist)
            max_index = hehelist.index(max_val)
            if max_index == 0:
                consensus += 'A'
            elif max_index == 1:
                consensus += 'C'
            elif max_index == 2:
                consensus += 'G'
            else:
                consensus += 'T'
        print(inFile + ': ')
        print(w_score)
        print(consensus + '\n')

    else:
        raise Usage("somthing went rly wrong sos go get help")
    
       
if __name__ == "__main__":
    main(inFile='p1860Crisprs', options=['-i 1000', '-k 13', '-p 1'])
    

