import networkx as nx
import matplotlib.pyplot as plt
import rdflib

filename = 'cleaned_maingraph.nt'
g = rdflib.Graph()
g.parse(filename, format="nt")

query = """
    PREFIX PMC: <https://pubmed.ncbi.nlm.nih.gov/>
    PREFIX EDAM: <http://edamontology.org/>
    PREFIX RDF: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX DCT: <http://purl.org/dc/terms/>
    PREFIX BIOLINK: <https://w3id.org/biolink/vocab/>
    PREFIX ENSEMBL: <http://identifiers.org/ensembl/>
    PREFIX NCBIGENE: <http://identifiers.org/ncbigene/>
                
    SELECT ?pmid ?gene ?logfold
    WHERE {
        ?pmid (<>|!<>)* ?gene .
        ?s BIOLINK:symbol | ENSEMBL:id | NCBIGENE:id ?gene . 
        ?pmid (<>|!<>)* ?logfold .
        ?s EDAM:data_3754 ?logfold . 
        FILTER(?pmid IN (PMC:28184278, PMC:32460837, PMC:33262327))
    }
"""

results = g.query(query)

# Create graph
subG = nx.Graph()
for result in results["results"]["bindings"]:
    subject = result["subject"]["value"]
    predicate = result["predicate"]["value"]
    object = result["object"]["value"]
    subG.add_edge(subject, object, label=predicate)

# Draw graph
pos = nx.spring_layout(subG)
nx.draw(subG, pos, with_labels=True, node_size=500, node_color="lightblue", font_size=8)
edge_labels = nx.get_edge_attributes(subG, 'label')
nx.draw_networkx_edge_labels(subG, pos, edge_labels=edge_labels)

plt.axis('off')
plt.tight_layout()
plt.show()





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
plt.plot()
#nx.draw(sg, with_labels=True, font_weight='bold')
#nx.draw(sg, node_size=1, node_color='blue', with_labels=False)
"""