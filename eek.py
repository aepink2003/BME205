import sys

class node:
    """ Represents a node in a Directed Acyclic Graph (DAG). """

    def __init__(self, label):
        """
        Initialize a new node with a label.
        
        Parameters:
        label (str): The label of the node.
        """
        # List of parent nodes and the weights of their edges.
        self.label = label
        self.parent_nodes = []  
        # List of target nodes and the weights of their edges.
        self.target_nodes = [] 
        # Boolean flag to mark if the node has been visited. 
        self.visited = False    

class DAG:
    """ Represents a Directed Acyclic Graph (DAG) and operations to perform on it. """

    def __init__(self):
        """
        Initialize a new DAG.
        """
        # Dictionary to store nodes with their labels as keys.
        self.nodes_dict = {} 
        # Dictionary to store the longest distance from source to each node. 
        self.distance = {} 
         # Dictionary to store the previous node for the longest path.   
        self.backtrack = {}  
    def add_node(self, label):
        """
        Add a new node to the DAG, or return it if it already exists.

        Parameters:
        label (str): The label of the node to add or retrieve.

        Returns:
        node: The new or existing node.
        """
        if label in self.nodes_dict:
            return self.nodes_dict[label]

        new_node = node(label)
        self.nodes_dict[label] = new_node
        return new_node

    def contruct_dag(self, adj_list_text):
        """
        Construct the DAG using an adjacency list text representation.

        Parameters:
        adj_list_text (list of str): The adjacency list where each element represents an edge.
        """
        for line in adj_list_text:
            nodeA, tmp = line.split("->")
            nodeB, weight = tmp.split(":")
            weight = int(weight)

            from_node = self.add_node(nodeA)
            to_node = self.add_node(nodeB)

            from_node.target_nodes.append((to_node, weight))
            to_node.parent_nodes.append((from_node, weight))

    def topo_sort_u(self, node, stack):
        """
        Utility function for topological sort.

        Parameters:
        node (node): The current node to process.
        stack (list of str): The stack to keep track of the sorted nodes.
        """
        node.visited = True
        for node2,_ in node.target_nodes:
            if not node2.visited:
                self.topo_sort_u(node2, stack)
        stack.insert(0, node.label)
        
    def topo_sort(self):
        """
        Perform topological sort on the DAG.

        Returns:
        list of str: The list of node labels in topologically sorted order.
        """
        stack = []
        for node in self.nodes_dict.values():
            if not node.visited:
                self.topo_sort_u(node, stack)
        return stack

    def long_path(self, source, sink):
        """
        Find the longest path in the DAG from a source node to a sink node.

        Parameters:
        source (str): The label of the source node.
        sink (str): The label of the sink node.

        Returns:
        tuple: A tuple containing the longest distance and the path as a list of node labels.
        """
        for label in self.nodes_dict:
            self.distance[label] = -float("Inf")

        self.distance[source] = 0
        self.backtrack[source] = None

        top_order = self.topo_sort()
        for label in top_order:
            current_node = self.nodes_dict[label]
            for v, weight in current_node.target_nodes:
                if self.distance[v.label] < self.distance[label] + weight:
                    self.distance[v.label] = self.distance[label] + weight
                    self.backtrack[v.label] = label

        path = [sink]
        curr = self.backtrack[sink]
        while curr != source:
            path = [curr] + path
            curr = self.backtrack[curr]
        path = [source] + path
        return self.distance[sink], path


if __name__ == "__main__":
    # Check if the input file name is provided.
    if len(sys.argv) < 2:
        print("Usage: python script.py filename")
        sys.exit(1)

    # The input file name.
    infile = sys.argv[1]

    # Open the input file to read the adjacency list text.
    with open(infile, 'r') as file:
        input_data = file.read().splitlines()

    # The source node label.
    source = input_data[0]
    # The sink node label.
    sink = input_data[1]
    # The adjacency list text.
    adj_list_text = input_data[2:]

    # Create and construct the DAG.
    graph = DAG()
    graph.contruct_dag(adj_list_text)

    # Find the longest path.
    longest_dist, long_path = graph.long_path(source, sink)

    # Output file name.
    outfile = infile + ".out"

    # Open the output file to write the results.
    with open(outfile, 'w') as file:
        print(longest_dist, file=file)
        print("->".join(long_path), file=file)

    print(f"Results written to {outfile}")
    # # Read the input from standard input.
    # input = sys.stdin.read().splitlines()
    # # The source node label.
    # source = input[0] 
    # # The sink node label. 
    # sink = input[1]    
    # adj_list_text = input[2:]  
    # # The adjacency list text.
    # graph = DAG()
    # graph.contruct_dag(adj_list_text)
    # longest_dist, long_path = graph.long_path(source, sink)
    # print(longest_dist)
    # print("->".join(long_path))
