def calculate_probability(x, pi, emission_matrix):
    probability = 1.0
    for i in range(len(x)):
        state = pi[i]
        symbol = x[i]
        probability *= emission_matrix[state][symbol]
    return probability

#  daata
x = "yzyxyzyzyzyyyzxzxxzyyyyzyxxxzzzyyzzzyyyyyzzyyxxxzx"
pi = "AABBBBABBABABBAAABBBABAAABBAABAABBAAABAABBABBBAABA"
emission_matrix = {
    'A': {'x': 0.13, 'y': 0.212, 'z': 0.659},
    'B': {'x': 0.181, 'y': 0.637, 'z': 0.183}
}

# Calculate the probability
probability = calculate_probability(x, pi, emission_matrix)
print(f"{probability}")