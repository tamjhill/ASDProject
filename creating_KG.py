import networkx as nx
from rdflib import Graph
import matplotlib.pyplot as plt

# Load the N-Triples file
g = Graph().parse("knowledge-graph.nt", format="nt")

# Create a NetworkX graph from the RDF graph
nx_graph = nx.DiGraph()
for subj, pred, obj in g:
    nx_graph.add_edge(str(subj), str(obj), label=str(pred))

# Draw the graph
pos = nx.spring_layout(nx_graph)
nx.draw(nx_graph, pos, with_labels=True)
edge_labels = nx.get_edge_attributes(nx_graph, 'label')
nx.draw_networkx_edge_labels(nx_graph, pos, edge_labels=edge_labels)

# Display the graph
plt.show()