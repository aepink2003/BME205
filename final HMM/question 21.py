import numpy as np

# Input Data
sequence = "zxzzzzxzxzxxxyxzyyzxxxzxyzxyyxyzzzyxxyyyyxzyxzzyxzzyzzxzzzyyyyxyxyyzxxyyxzyzyzxzxyzzxyxyzxzzzxzxxyxy"
states = ["A", "B"]
alphabet = {"x": 0, "y": 1, "z": 2}

# Transition Matrix
transition_matrix = np.array([
    [0.86, 0.14],
    [0.41, 0.59]
])

# Emission Matrix
emission_matrix = np.array([
    [0.669, 0.312, 0.02],
    [0.467, 0.295, 0.238]
])

# Initialize Viterbi Matrix
n_states = len(states)
n_sequence = len(sequence)
viterbi_matrix = np.zeros((n_states, n_sequence))
backtrack = np.zeros((n_states, n_sequence), dtype=int)

# Initialization step
viterbi_matrix[:, 0] = 1.0 * emission_matrix[:, alphabet[sequence[0]]]

# Viterbi algorithm
for i in range(1, n_sequence):
    for j in range(n_states):
        prob = viterbi_matrix[:, i-1] * transition_matrix[:, j] * emission_matrix[j, alphabet[sequence[i]]]
        viterbi_matrix[j, i] = np.max(prob)
        backtrack[j, i] = np.argmax(prob)

# Backtracking
path = np.zeros(n_sequence, dtype=int)
path[-1] = np.argmax(viterbi_matrix[:, -1])

for i in range(n_sequence-2, -1, -1):
    path[i] = backtrack[path[i+1], i+1]

# Converting numerical path to state symbols
decoded_path = "".join(states[state] for state in path)
print(decoded_path)
