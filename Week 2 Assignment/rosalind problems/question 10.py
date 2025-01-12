def overlap_graph(patterns):
    '''Initialize an empty dictionary to store the adjacency list.'''
    # Each key in this dictionary is a k-mer, and its value is a list of k-mers that it overlaps with.
    adjacency_list = {}
    
    # Initialize each k-mer's adjacency list as an empty list.
    for pattern in patterns:
        adjacency_list[pattern] = []

    # For each pair of k-mers (pattern and pattern_prime) in the collection:
    for pattern in patterns:
        for pattern_prime in patterns:
            # Check if the two k-mers are different and if the suffix of 'pattern' 
            # matches the prefix of 'pattern_prime'.
            if pattern != pattern_prime and pattern[1:] == pattern_prime[:-1]:
                # If they overlap, add 'pattern_prime' to the adjacency list of 'pattern'.
                adjacency_list[pattern].append(pattern_prime)
                
    return adjacency_list

def print_adjacency_list(adjacency_list):
    # For each k-mer and its overlaps in the adjacency list:
    for key, values in adjacency_list.items():
        # For each overlapping k-mer:
        for value in values:
            # Print the k-mer and its overlapping k-mer in the desired format.
            print(f"{key} -> {value}")

# Sample Dataset
with open('rosalind_ba3c.txt', 'r') as file:
    patterns = file.read().splitlines()


# Generate the overlap graph using the 'overlap_graph' function.
adj_list = overlap_graph(patterns)

# Print the overlap graph using the 'print_adjacency_list' function.
print_adjacency_list(adj_list)
