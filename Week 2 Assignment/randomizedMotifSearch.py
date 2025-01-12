import numpy as np
import math
import sys

class FastAreader:
    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''

        self.fname = fname
        self.fileH = None

    def doOpen(self):
        if self.fname == None:
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
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
        )

        self.parser.add_argument('-i', '--iterations', nargs='?', type=int, default=1000, action='store',
                                 help='number of iterations')
        self.parser.add_argument('-k', '--motifLength', nargs='?', type=int, default=13, action='store',
                                 help='length of motif')
        self.parser.add_argument('-p', '--pseudoCount', nargs='?', type=float, default=1, action='store',
                                 help='pseudocount to add to each nucleotide count')
        self.parser.add_argument('-m', '--showMotifs', action='store_true',
                                 help='Show motifs in FastA Sequence in output')
        self.parser.add_argument('-g', '--gibbsSample', action='store_true',
                                 help='Gibbs Sampling instead of random sampling')
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)
class Motif:
    """
    
    """
    def __init__(self, sequences, psCount = 1,motifLen = 13):
            '''Construct an object for finding the best motifs from a set of given sequences.'''
            # Contains a list of sequences
            self.seqs = sequences                        
            self.k = int(motifLen)     
            # Holds pseudocount value         
            self.pseudo = float(psCount)       
            allSeqs = ''.join(self.seqs)
            nucleotides = ['A', 'C', 'T', 'G']
            self.nucleotide_counts = {nt: (self.pseudocount + allSeqs.count(nt)) for nt in nucleotides}
            self.bestMotifs = []        
            self.bestProfile = {}
            # Initialize with lowest possible score
            self.bestScore = 0   
    def generateRandomMotifs(self, sequences):
        '''Generate and return random motifs from the provided sequences.'''
        # mpty list for storing motifs
        motifs = []                                            
        for sequence in sequences:                           
            startIdx = np.random.randint(0, len(sequence) - self.motifLength + 1) 
            motifs.append(sequence[startIdx:startIdx + self.motifLength])  
        return motifs
    def makeProfile(self, motifs):
        '''Generate and return a nucleotide count profile from the provided motifs.'''
        # init a profile for each position in the motif
        profile = {position: {} for position in range(self.motifLength)}  
        for position in profile:           
            # Extract nucleotides at the current position from all motifs                             
            nucleotides = [motif[position] for motif in motifs]         
            profile[position] = {nucleo: (self.smoothing + nucleotides.count(nucleo)) for nucleo in ['A', 'C', 'T', 'G']}  
        return profile
    def computeRelativeEntropy(self, profileMatrix):
        '''Compute and return the relative entropy score of a nucleotide profile matrix.'''
        relativeEntropy = 0      
        # Iterate over each position in the motif                                 
        for position in range(self.motifLength):                  
            for nucleotide in ['A', 'C', 'T', 'G']:
                freq = profileMatrix[position][nucleotide] / sum(profileMatrix[position].values())  
                baseProb = self.baseFrequency[nucleotide] / sum(self.baseFrequency.values())      
                # Increment the relative entropy score
                relativeEntropy += freq * math.log2(freq/baseProb)  
        return relativeEntropy
    def BestMotifs(self, sequences, profileMatrix):
        '''Determine and return motifs from sequences that best match a provided nucleotide profile matrix.'''
        # Start with an empty list for the selected motifs
        selectedMotifs = []                                       
        for sequence in sequences:
            highestProbability = 0                                 
            bestMatchingMotif = ''
            
            #need help here idk
            for index in range(len(sequence) - self.motifLength + 1):  
                currentMotif = sequence[index:index+self.motifLength]
                 # Initialize joint probability for the current motif
                motifProbability = 1                              
               
                    
                
            # Add the best matching motif for the current sequence to the list               
        return selectedMotifs
    
    def identifyOptimalMotif(self, maxIterations=1000, useGibbs=False):
    '''Optimize the motif selection, profile, and score from a set of sequences.
    This method uses an Expectation Maximization Algorithm through multiple iterations to identify the
    optimal motif. The user has the option to utilize either Gibbs sampling or random sampling.'''
    
    for iteration in range(maxIterations):
        if useGibbs and iteration != 0:  # If Gibbs sampling is selected, start with random sampling for the first iteration
            currentMotifs = self.sampleMotifWithGibbs(prevMotifs)  # Generate motifs based on the previous set
            prevMotifs = currentMotifs  # Store motifs for the next Gibbs sampling round
            currentProfile = self.generateProfileMatrix(currentMotifs)  # Generate a profile from the current motifs
            currentScore = self.calculateProfileScore(currentProfile)   # Compute the score for the current profile
        else:
            currentMotifs = self.generateRandomMotifs(self.sequences)  # Generate a random set of motifs
            prevMotifs = currentMotifs  # Store the initial set for Gibbs sampling
            currentProfile = self.generateProfileMatrix(currentMotifs)  # Generate profile from the current motifs
            currentScore = self.calculateProfileScore(currentProfile)   # Compute the score for the current profile
        
        iterationBestMotifs = currentMotifs  # Store the motifs, profile, and score to identify the best in each iteration
        iterationBestProfile = currentProfile
        iterationBestScore = currentScore
        
        while True:
            currentMotifs = self.selectMotifs(self.sequences, currentProfile)  # Get a new set of motifs based on the current profile
            currentProfile = self.generateProfileMatrix(currentMotifs)  # Generate a profile from the current motifs
            currentScore = self.calculateProfileScore(currentProfile)   # Compute the score for the current profile
            
            if currentScore > iterationBestScore:  # If the current score is better than the iteration's best, update the iteration's best values
                iterationBestMotifs = currentMotifs
                iterationBestProfile = currentProfile
                iterationBestScore = currentScore
            else:  # If the score does not improve or gets worse, we've found the best motifs for this iteration
                if iterationBestScore > self.overallBestScore:  # If the iteration's best score is better than the overall best, update the overall best values
                    self.overallBestMotifs = iterationBestMotifs
                    self.overallBestProfile = iterationBestProfile
                    self.overallBestScore = iterationBestScore
                break  # End the current iteration and proceed to the next one

    
    
    
    
    
    
    


 
