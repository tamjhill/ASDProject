import rdflib
from rdflib import Graph, URIRef, Literal, BNode
import requests
from rdflib.plugins.sparql import prepareQuery
import networkx as nx
import matplotlib.pyplot as plt
from rdflib.extras.external_graph_libs import rdflib_to_networkx_multidigraph

filename = "test_graph4.ttl"

g = rdflib.Graph()
g.parse(filename)
print(len(g))
query = prepareQuery("""
    SELECT DISTINCT ?doi
    WHERE {
        ?doi (<>|!<>)* ?object .
        FILTER(CONTAINS(STR(?object), "CHN1")) .
        FILTER(REGEX(STR(?doi), "^https://doi.org/"))             
        }
""")
results = g.query(query)
for row in results:
    print(row)
print(len(results))

#?intermediate ?predicate ?object .
"""

"""