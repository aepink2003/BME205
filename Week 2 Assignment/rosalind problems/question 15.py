import sys


class node:
    def __init__(self, label):
        self.label = label
        self.parent_nodes = []
        self.target_nodes = []
        self.visited = False


class DAG:
    def __init__(self):
        self.nodes_dict = {}
        self.distance = {}
        self.backtrack = {}

    def add_node(self, label):
        if label in self.nodes_dict:
            return self.nodes_dict[label]

        new_node = node(label)
        self.nodes_dict[label] = new_node
        return new_node

    def contruct_dag(self, adj_list_text):
        for line in adj_list_text:
            nodeA, tmp = line.split("->")
            nodeB, weight = tmp.split(":")
            weight = int(weight)

            from_node = self.add_node(nodeA)
            to_node = self.add_node(nodeB)

            from_node.target_nodes.append((to_node, weight))
            to_node.parent_nodes.append((from_node, weight))

    def topological_sort_u(self, node, stack):
        node.visited = True
        for node2,_ in node.target_nodes:
            if not node2.visited:
                self.topological_sort_u(node2, stack)
        stack.insert(0, node.label)

    def topological_sort(self):
        stack = []
        for node in self.nodes_dict.values():
            if not node.visited:
                self.topological_sort_u(node, stack)
        return stack

    def long_path(self, source, sink):
        for label in self.nodes_dict:
            self.distance[label] = -float("Inf")

        self.distance[source] = 0
        self.backtrack[source] = None

        top_order = self.topological_sort()
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
    
    input = sys.stdin.read().splitlines()
    source = input[0]
    sink = input[1]
    adj_list_text = input[2:]

    graph = DAG()
    graph.contruct_dag(adj_list_text)
    longest_dist, long_path = graph.long_path(source, sink)
    print(longest_dist)
    print("->".join(long_path))