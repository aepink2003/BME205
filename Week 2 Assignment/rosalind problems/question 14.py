# def string_reconstruction(k, patterns):
#     """Reconstructs a string from a set of k-mers.

#     Args:
#         k: The length of the k-mers.
#         patterns: A list of k-mers.

#     Returns:
#         A string reconstructed from the k-mers.
#     """

#     # Create the adjacency list
#     adj_list = {}
#     for pattern in patterns:
#         prefix = pattern[:-1]
#         suffix = pattern[1:]
#         if prefix in adj_list:
#             adj_list[prefix].append(suffix)
#         else:
#             adj_list[prefix] = [suffix]

#     # Find the start and end nodes
#     all_nodes = set(adj_list.keys()).union(set([x for v in adj_list.values() for x in v]))
#     in_degrees = {node: 0 for node in all_nodes}
#     out_degrees = {node: 0 for node in all_nodes}
#     for node in all_nodes:
#         if node in adj_list:
#             out_degrees[node] = len(adj_list[node])
#             for neighbor in adj_list[node]:
#                 in_degrees[neighbor] += 1
#         else:
#             out_degrees[node] = 0
#     start_node = None
#     for node in all_nodes:
#         if in_degrees[node] < out_degrees[node]:
#             start_node = node
#             break
#     if not start_node:
#         start_node = list(all_nodes)[0]

#     # Follow the Eulerian path
#     path = [start_node]
#     while len(adj_list) > 0:
#         current_node = path[-1]
#         if current_node in adj_list:
#             next_node = adj_list[current_node][0]
#             path.append(next_node)
#             if len(adj_list[current_node]) == 1:
#                 del adj_list[current_node]
#             else:
#                 adj_list[current_node] = adj_list[current_node][1:]
#         else:
#             break

#     # Construct the string from the path
#     string = path[0]
#     for i in range(1, len(path)):
#         string += path[i][-1]

#     # Return the string without the first and last k-mer
#     return string[k-1:]

# # Example usage:

# # k = 4
# # patterns = []
# # with open("rosalind_ba3h-4.txt") as f:
# #     for line in f:
# #         patterns.append(line.strip())

# # # Add a separator between the input and output
# # output = "///\n" + string_reconstruction(k, patterns)
# # print(output)

# def something_something(filename):
#     graph = {}
#     with open(filename, 'r') as f:
#         for line in f:
#             parts = line.strip().split(' -> ')
#             key = int(parts[0])
#             values = list(map(int, parts[1].split(',')))
#             graph[key] = values
#     return graph

from collections import defaultdict

def de_bruijn_graph(k, patterns):
    graph = defaultdict(list)
    
    for pattern in patterns:
        prefix = pattern[:k-1]
        suffix = pattern[1:]
        graph[prefix].append(suffix)
    
    return graph

def eulerian_path(graph):
    def find_starting_node(graph):
        in_degree = defaultdict(int)
        out_degree = defaultdict(int)

        for node, neighbors in graph.items():
            out_degree[node] = len(neighbors)
            for neighbor in neighbors:
                in_degree[neighbor] += 1

        start_node = end_node = None

        for node in out_degree.keys():
            if out_degree[node] > in_degree[node]:
                start_node = node
            elif out_degree[node] < in_degree[node]:
                end_node = node

        return start_node, end_node

    def find_eulerian_path(node, graph):
        path = [node]

        while node in graph:
            if len(graph[node]) > 0:
                next_node = graph[node].pop()
                path.append(next_node)
                node = next_node
            else:
                break

        return path

    start_node, end_node = find_starting_node(graph)
    graph[end_node].append(start_node)

    path = find_eulerian_path(start_node, graph)

    return path

def reconstruct_string_from_kmers(k, patterns):
    graph = de_bruijn_graph(k, patterns)
    path = eulerian_path(graph)
    reconstructed_string = path[0] + ''.join(node[-1] for node in path[1:])
    
    return reconstructed_string

# Reads the input from the 'dataset.txt' file
with open('rosalind_ba3h-4.txt', 'r') as f:
    k = int(f.readline().strip())
    patterns = [line.strip() for line in f.readlines()]

# Calling reconstruction function
result = reconstruct_string_from_kmers(k, patterns)

# Print the result
print(result)

