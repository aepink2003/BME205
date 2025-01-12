import numpy as np
import sys

class HMM:
    def __init__(self, states, transition, emission):
        self.states = states
        self.transition = transition
        self.emission = emission

    def viterbi(self, sequence, alphabet):
        # Initialize matrices
        T = len(sequence)
        S = len(self.states)
        viterbi_matrix = np.zeros((S, T))
        path_matrix = np.zeros((S, T), dtype=int)

        # Initialize first column
        for s in range(S):
            viterbi_matrix[s][0] = self.emission[self.states[s]][sequence[0]]

        # Iterate over the sequence
        for t in range(1, T):
            for s in range(S):
                (prob, state) = max((viterbi_matrix[s_prime][t-1] *
                                     self.transition[self.states[s_prime]][self.states[s]] *
                                     self.emission[self.states[s]][sequence[t]], s_prime) 
                                    for s_prime in range(S))
                viterbi_matrix[s][t] = prob
                path_matrix[s][t] = state

        # Backtrack to find the path
        best_path = [np.argmax(viterbi_matrix[:, T-1])]
        for t in range(T-1, 0, -1):
            best_path.insert(0, path_matrix[best_path[0]][t])

        return ''.join(self.states[i] for i in best_path)

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
    path = hmm.viterbi(sequence, alphabet)
    
    with open(inFile + '.out', 'w') as f:
        print(path, file=f)

if __name__ == "__main__":
    main(inFile='rosalind_ba10c.txt') 
