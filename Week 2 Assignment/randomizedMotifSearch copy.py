import argparse
import random
import math

class RandomizedMotifFinder:
    """
    A class to represent a randomized motif search on sequences.
    """

    def __init__(self, sequences, motifLength, maxIterations, pseudoCount):
        """
        Initializes the RandomizedMotifFinder with sequences, motif length, iterations, and pseudocount.
        """
        self.sequences = sequences
        self.motifLength = motifLength
        self.maxIterations = maxIterations
        self.pseudoCount = pseudoCount
        self.optimalMotifs = []
        self.optimalScore = float('inf')

    def _initializeMotifsRandomly(self):
        """
        Initialize the motifs from the sequences randomly.
        """
        motifs = []
        for seq in self.sequences:
            startPos = random.randint(0, len(seq) - self.motifLength)
            motifs.append(seq[startPos:startPos+self.motifLength])
        return motifs

    def _generateProfileMatrix(self, motifs):
        """
        Generates a profile matrix (frequencies) from given motifs.
        """
        # Initialize the profile with pseudocounts
        profile = {'A': [self.pseudoCount]*self.motifLength, 
                   'C': [self.pseudoCount]*self.motifLength,
                   'G': [self.pseudoCount]*self.motifLength, 
                   'T': [self.pseudoCount]*self.motifLength}
        
        # Increment the counts based on motifs
        for motif in motifs:
            for idx, nucleotide in enumerate(motif):
                profile[nucleotide][idx] += 1
        
        # Convert counts to frequencies
        totalMotifs = len(motifs) + 4 * self.pseudoCount
        for key, values in profile.items():
            profile[key] = [val / totalMotifs for val in values]
        return profile

    def _selectProbableMotif(self, sequence, profile):
        """
        Select the most probable motif from a sequence based on a profile.
        """
        maxProbability = -1
        bestMotif = sequence[:self.motifLength]
        
        # Check each possible motif in the sequence
        for i in range(len(sequence) - self.motifLength + 1):
            subSeq = sequence[i:i+self.motifLength]
            probability = 1
            for j, nucleotide in enumerate(subSeq):
                probability *= profile[nucleotide][j]
            if probability > maxProbability:
                maxProbability = probability
                bestMotif = subSeq
        return bestMotif

    def _calculateMotifScore(self, motifs, profile):
        """
        Calculate the score of motifs based on the profile.
        """
        score = 0
        for i in range(self.motifLength):
            column = [profile[nucleotide][i] for nucleotide in "ACTG"]
            maxFrequency = max(column)
            score += (1 - maxFrequency)
        return score

    def findMotif(self):
        """
        Perform the randomized motif search.
        """
        optimalMotifs = self._initializeMotifsRandomly()
        optimalScore = self._calculateMotifScore(optimalMotifs, self._generateProfileMatrix(optimalMotifs))
        
        # Iterate to improve the motifs
        for _ in range(self.maxIterations):
            motifs = self._initializeMotifsRandomly()
            while True:
                currentProfile = self._generateProfileMatrix(motifs)
                newMotifs = [self._selectProbableMotif(seq, currentProfile) for seq in self.sequences]
                score = self._calculateMotifScore(newMotifs, currentProfile)
                
                # If find better motifs, update our best ones
                if score < optimalScore:
                    optimalScore = score
                    optimalMotifs = newMotifs[:]
                else:
                    break
            
            # Update optimal score and motifs
            if optimalScore < self.optimalScore:
                self.optimalScore = optimalScore
                self.optimalMotifs = optimalMotifs[:]
        return self.optimalMotifs, self.optimalScore

    def consensusSequence(self, motifs):
        """
        Derive a consensus sequence from motifs.
        """
        consensus = ''
        for i in range(self.motifLength):
            column = [motif[i] for motif in motifs]
            consensus += max(set(column), key=column.count)
        return consensus

if __name__ == "__main__":
    # Argument parsing for command line usage
    parser = argparse.ArgumentParser(description="Determine consensus motif from sequences using randomized motif search.")
    parser.add_argument('-k', type=int, required=True, help="Length of motif.")
    parser.add_argument('-i', type=int, required=True, help="Number of iterations.")
    parser.add_argument('-p', type=float, default=0.1, help="Pseudocount value.")
    args = parser.parse_args()

    # Reading sequences
    sequences = []  
    try:
        while True:
            header = input().strip()
            sequence = input().strip()
            sequences.append(sequence)
    except EOFError:
        pass

    # Perform the motif search and display results
    motifFinder = RandomizedMotifFinder(sequences, args.k, args.i, args.p)
    bestMotifs, bestScore = motifFinder.findMotif()
    print(f"Best motifs: {bestMotifs}")
    print(f"Consensus: {motifFinder.consensusSequence(bestMotifs)}")
    print(f"Score: {bestScore}")
