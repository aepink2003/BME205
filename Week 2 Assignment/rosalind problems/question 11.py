def de_bruijn_graph(k, text):
   edges = {}
   for i in range(len(text) - k + 1):
       kmer = text[i:i + k]
       prefix = kmer[:k - 1]
       suffix = kmer[1:]
       if prefix in edges:
           edges[prefix].append(suffix)
       else:
           edges[prefix] = [suffix]
   return edges

def format_output(edges):
    sorted_keys = sorted(edges.keys())
    answer = []
    for key in sorted_keys:
        answer.append(f"{key} -> {','.join(edges[key])}")
    return "\n".join(answer)

def input(filename):
    with open(filename, 'r') as f:
        k = int(f.readline().strip())
        text = f.readline().strip()
    return k, text


# filename = "rosalind_ba3d-13.txt"
# k, text = input(filename)
# edges = de_bruijn_graph(k, text)
# print(format_output(edges))

with open('output.txt','w') as x:
    filename = "rosalind_ba3d-13.txt"
    k, text = input(filename)
    edges = de_bruijn_graph(k, text)
    x.write(format_output(edges))
 
# # Allowing user input for k and text
# k = int(input("Enter the value of k: "))
# text = input("Enter the string Text: ")
 
# # Construct the de Bruijn graph
# graph = de_bruijn_graph(k, text)
 
# # Convert the graph into the required format
# for key in graph:
#    suffixes = ",".join(graph[key])
#    print(f"{key} -> {suffixes}")
