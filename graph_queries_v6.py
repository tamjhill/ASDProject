from cycler import V
import rdflib
from rdflib import Graph, URIRef, Literal, BNode
import requests
from rdflib.plugins.sparql import prepareQuery
import networkx as nx
import matplotlib.pyplot as plt
from rdflib.extras.external_graph_libs import rdflib_to_networkx_multidigraph

filename = "test_graph2.nt"

g = rdflib.Graph()
g.parse(filename)
#print(len(g))
test_query = ("""
    PREFIX EDAM: <http://edamontology.org/>
    PREFIX RDF: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX DCT: <http://purl.org/dc/terms/>
    PREFIX BIOLINK: <https://w3id.org/biolink/vocab/>
    
    SELECT DISTINCT ?pmid ?gene
    WHERE {
        ?pmid (<>|!<>)* ?gene .
        ?pmid RDF:type DCT:identifier .
        FILTER(CONTAINS(STR(?gene), "SOX4")) .             
        }
""")
results = g.query(test_query)
for row in results:
    print(f"Article {row.pmid} includes gene {row.gene}")

#?intermediate ?predicate ?object .
"""

"""