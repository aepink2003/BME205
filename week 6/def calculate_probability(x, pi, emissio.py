def calculate_probability(x, pi, emission_matrix):
    probability = 1.0
    for i in range(len(x)):
        state = pi[i]
        symbol = x[i]
        probability *= emission_matrix[state][symbol]
    return probability

# Sample dataset
x = "xxyzyxzzxzxyxyyzxxzzxxyyxxyxyzzxxyzyzxzxxyxyyzxxzx"
pi = "BBBAAABABABBBBBBAAAAAABAAAABABABBBBBABAABABABABBBB"
emission_matrix = {
    'A': {'x': 0.612, 'y': 0.314, 'z': 0.074},
    'B': {'x': 0.346, 'y': 0.317, 'z': 0.336}
}

# Calculate the probability
probability = calculate_probability(x, pi, emission_matrix)
print(f"The conditional probability Pr(x|π) is: {probability:.2e}")