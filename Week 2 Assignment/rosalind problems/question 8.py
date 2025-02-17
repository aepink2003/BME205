def k_mer_composition(k, text):
    # Create an empty list to store the k-mers
    k_mers = []

    # Determine the number of k-mers by subtracting k from the length of the text and adding 1
    num_k_mers = len(text) - k + 1

    # Iterate through the text to extract all k-mers
    for i in range(num_k_mers):
        # Extract a k-mer from the current position to the current position + k
        k_mer = text[i:i+k]

        # Append the k-mer to the list of k-mers
        k_mers.append(k_mer)

    # Sort the list of k-mers (this step is optional and is done here to match the example provided)
    k_mers.sort()

    return k_mers

# Sample data
k = 5
text = 'CAATCCAAC'

# Get the k-mer composition
composition_k = k_mer_composition(k, text)

# Output the k-mer composition
for k_mer in composition_k:
    print(k_mer)
