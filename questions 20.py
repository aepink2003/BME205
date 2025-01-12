import numpy as np
import sys

def calculate_conditional_probability(sequence, hidden_path, states, emission_matrix):
    probability = 1.0
    for i, symbol in enumerate(sequence):
        state = hidden_path[i]
        probability *= emission_matrix[state][symbol]
    return probability

def main(inFile=None):
    with open(inFile, 'r') as fh:
        sequence = fh.readline().strip()
        fh.readline()  # Skip the delimiter
        alphabet = fh.readline().strip().split()
        hidden_path = fh.readline().strip()
        fh.readline()  # Skip the delimiter
        states = fh.readline().strip().split()
        fh.readline()  # Skip the delimiter
        emission_matrix = {state: {} for state in states}
        for state in states:
            values = list(map(float, fh.readline().strip().split()[1:]))
            emission_matrix[state] = dict(zip(alphabet, values))

    prob = calculate_conditional_probability(sequence, hidden_path, states, emission_matrix)
    
    with open(inFile + '.out', 'w') as f:
        print(prob, file=f)

if __name__ == "__main__":
    main(inFile='rosalind_ba10b.txt') 
