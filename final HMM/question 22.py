import numpy as np
import sys

class HMM:
    def __init__(self, states, transition, emission):
        self.states = states
        self.transition = transition
        self.emission = emission

    def probability_of_emission(self, sequence, alphabet):
        # Initialize the probabilities matrix
        prob_matrix = np.zeros((len(self.states), len(sequence)))

        # Initialize first column of the matrix
        for i, state in enumerate(self.states):
            prob_matrix[i][0] = self.emission[state][sequence[0]]

        # Fill in the matrix
        for j in range(1, len(sequence)):
            for i, state in enumerate(self.states):
                prob_matrix[i][j] = sum(prob_matrix[k][j-1] * self.transition[self.states[k]][state] * self.emission[state][sequence[j]] for k in range(len(self.states)))

        # Sum the probabilities for the last character
        return sum(prob_matrix[i][-1] for i in range(len(self.states)))

def main(inFile=None):
    with open(inFile, 'r') as fh:
        sequence = fh.readline().strip()
        fh.readline()  # Skip the delimiter
        alphabet = fh.readline().strip().split()
        fh.readline()  # Skip the delimiter
        states = fh.readline().strip().split()
        fh.readline()  # Skip the delimiter
        transition_matrix = {state: {} for state in states}
        for state in states:
            values = list(map(float, fh.readline().strip().split()[1:]))
            transition_matrix[state] = dict(zip(states, values))
        fh.readline()  # Skip the delimiter
        emission_matrix = {state: {} for state in states}
        for state in states:
            values = list(map(float, fh.readline().strip().split()[1:]))
            emission_matrix[state] = dict(zip(alphabet, values))

    hmm = HMM(states, transition_matrix, emission_matrix)
    prob = hmm.probability_of_emission(sequence, alphabet)
    
    with open(inFile + '.out', 'w') as f:
        print(prob, file=f)

if __name__ == "__main__":
    main(inFile='rosalind_ba10d.txt')
