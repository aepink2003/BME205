def string_spelled_by_genome_path(kmers):
    # Step 1: Initialize an empty string genome
    genome = kmers[0]

    # Step 3: For each remaining k-mer, add the last symbol of the k-mer to genome
    for kmer in kmers[1:]:
        genome += kmer[-1]

    # Step 4: Return genome
    return genome

# Sample dataset
with open('rosalind_ba3b-2.txt', 'r') as file:
    kmers = file.read().splitlines()

# Get the string spelled by the genome path
spelled_string = string_spelled_by_genome_path(kmers)

# Output the spelled string
print(spelled_string)
