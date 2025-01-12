def de_bruijn_graph(patterns):
   # Create an empty dictionary to represent the adjacency list
   graph = {}
   # Iterate through each k-mer in the collection of patterns
   for pattern in patterns:
       k = len(pattern)
       prefix = pattern[:-1]
       suffix = pattern[1:]
       # If the prefix is already in the graph, add the suffix to its list of neighbors
       if prefix in graph:
           graph[prefix].append(suffix)
       else:
           # If the prefix is not in the graph, create a new entry with the suffix as the neighbor
           graph[prefix] = [suffix]
   # Create the de Bruijn graph in the form of an adjacency list
   adjacency_list = []
   for key, values in graph.items():
       adjacency_list.append(f"{key} -> {', '.join(values)}")
   return adjacency_list
# Sample dataset

with open('rosalind_ba3e.txt', 'r') as f:
    patterns = [line.strip() for line in f.readlines()]
    
# patterns = [
    
#    "GAGG",
#    "CAGG",
#    "GGGG",
#    "GGGA",
#    "CAGG",
#    "AGGG",
#    "GGAG",
   
# ]
result = de_bruijn_graph(patterns)
for line in result:
   print(line)
 