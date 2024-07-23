




"""from fileinput import filename
from rdflib.extras.external_graph_libs import *
from rdflib import Graph, URIRef, Literal
import networkx as nx

graph = Graph()
filename = 
graph.parse(filename, format=“nt”)
nx_graph = rdflib_to_networkx_multidigraph(graph)
nx.write_graphml(nx_graph)


G = nx.MultiDiGraph()
for subj, obj, rel in triples_list:
    G.add_edge(subj, obj, label=rel) # subj = node 1, obj = node 2, rel = relation/edge name

pos = nx.spring_layout(G)

pos=pos
"""