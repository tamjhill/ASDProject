import rdflib
from rdflib import Graph, URIRef, Literal, BNode
import requests
import re
from rdflib.plugins.sparql import prepareQuery
import networkx as nx
import matplotlib.pyplot as plt
from rdflib.extras.external_graph_libs import rdflib_to_networkx_multidigraph


filename = "test_graph4.ttl"

g = rdflib.Graph()
g.parse(filename)
#print(len(g))
query = """
SELECT ?subject ?pvalue ?value
WHERE {
  ?subject ?pvalue ?value .
}
"""

results = g.query(query)

pvalue_pattern = re.compile(r'p[-\s]?value|pvalue', re.IGNORECASE)

def find_root(g, node, visited=None):
    if visited is None:
        visited = set()
    visited.add(node)
    parents = list(g.subjects(predicate=None, object=node))
    if not parents:
        return node
    for parent in parents:
        if parent not in visited:
            return find_root(g, parent, visited)
    return node

for row in results:
    predicate_str = str(row.pvalue)
    if pvalue_pattern.search(predicate_str):
        try:
            value = float(row.value)
            if value < 0.05:
                subject = row.subject
                root = find_root(g, subject)
                print(f"Root: {root}")
                print(f"Subject: {subject}")
                print(f"P-value predicate: {row.pvalue}")
                print(f"P-value: {value}")
                print()
        except ValueError:
            pass