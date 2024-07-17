import rdflib
from rdflib import Graph, URIRef, Literal, BNode
import requests
from rdflib.plugins.sparql import prepareQuery
import networkx as nx
import matplotlib.pyplot as plt
from rdflib.extras.external_graph_libs import rdflib_to_networkx_multidigraph

g = rdflib.Graph()
g.parse(filename)
print(len(g))
query = ("""
    PREFIX EDAM: <http://edamontology.org/>
    PREFIX RDF: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX DCT: <http://purl.org/dc/terms/>
    PREFIX BIOLINK: <https://w3id.org/biolink/vocab/>

    SELECT ?pmid ?gene ?pvalue
    WHERE {
        ?pmid (<>|!<>) ?gene .
        ?pmid (<>|!<>) ?pvalue .
        FILTER(REGEX(STR(?pvalue), "p[-\\s]?value|pvalue"))   .
        FILTER (?value < 0.05)          
    }
    ORDER BY ASC(?pvalue)
""")
results = g.query(query)
for row in results:
    print(row)
print(len(results))

#    GROUP BY ?pmid
#?intermediate ?predicate ?object .
"""

"""

filename = "test_graph1.ttl"

g = rdflib.Graph()
g.parse(filename)
query = """
    SELECT ?s ?p ?o
    WHERE {
        ?s ?p ?o
    }
    LIMIT 5
"""
results = g.query(query)
for row in results:
    print(row)