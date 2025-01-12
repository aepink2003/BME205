import numpy as np
import sys

def calculate_path_probability(hidden_path, states, transition_matrix):
    # Assuming equal initial probabilities
    initial_prob = 1.0 / len(states)
    probability = initial_prob

    for i in range(1, len(hidden_path)):
        prev_state = hidden_path[i - 1]
        current_state = hidden_path[i]
        probability *= transition_matrix[prev_state][current_state]

    return probability

def main(inFile=None):
    with open(inFile, 'r') as fh:
        hidden_path = fh.readline().strip()
        fh.readline()  # Skip the delimiter
        states = fh.readline().strip().split()
        fh.readline()  # Skip the delimiter
        transition_matrix = {state: {} for state in states}
        for state in states:
            values = list(map(float, fh.readline().strip().split()[1:]))
            transition_matrix[state] = dict(zip(states, values))

    prob = calculate_path_probability(hidden_path, states, transition_matrix)
    
    with open(inFile + '.out', 'w') as f:
        print(prob, file=f)

if __name__ == "__main__":
    main(inFile='rosalind_ba10a.txt') 
