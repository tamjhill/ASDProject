"""ignore for now - but could use serialize method for combing with PrimeKG"""


import networkx as nx
import rdflib
import matplotlib.pyplot as plt
from rdflib.extras.external_graph_libs import rdflib_to_networkx_multidigraph
import streamlit as st
import networkx as nx
import plotly.graph_objects as go
import dash
import plotly.graph_objects as go

# Create a new RDFLib Graph adding each triple file
graph = rdflib.Graph()
graph.parse("result3.nt", format="nt")
graph.parse("result4.nt", format="nt")
graph.parse("result5.nt", format="nt")

# Optionally, serialize the combined graph to a new file
graph.serialize("combined_graphv2.nt", format="nt")

# SPARQL query
q = """
    PREFIX has_output: <http://edamontology.org/has_output>

    SELECT ?s ?p ?o
    WHERE {
        ?s ?p ?o .
        FILTER (CONTAINS(str(?o), "FUNCTION??"))
}
""" 

# Execute the SPARQL query on the RDFLib Graph
qres = graph.query(q)

# Iterate over the query results and print the matching triples
#for row in qres:
#    print(row)

#convert to nx graph for plotting
G = rdflib_to_networkx_multidigraph(qres)

# Plot nx instance of RDF Graph
pos = nx.spring_layout(G, scale=2)
edge_labels = nx.get_edge_attributes(G, 'r')
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
nx.draw(G, with_labels=True)
plt.show()


# Extract the nodes (subjects and objects) and edges (predicates)
#nodes = list(set([str(s) for s, _, _ in qres] + [str(o) for _, _, o in qres]))
#edges = []
#for s, p, o in qres:
#    edges.append((str(s), str(p), str(o)))

# Create the network graph
#edge_trace = []

#for edge in edges:
#    x0, y0 = nodes.index(edge[0]), nodes.index(edge[2])
#    x1, y1 = nodes.index(edge[2]), nodes.index(edge[0])
#    edge_x = [x0, x1, None]
#    edge_y = [y0, y1, None]
#    edge_trace.append(
#        go.Scatter(
#            x=edge_x,
#            y=edge_y,
#            mode="lines",
#            line=dict(color="rgb(210,210,210)", width=1),
#            hoverinfo="none"
#        )
#    )

#node_positions = [(nodes.index(node), nodes.index(node)) for node in nodes]

#node_trace = go.Scatter(
#    x=[pos[0] for pos in node_positions],
#    y=[pos[1] for pos in node_positions],
#    text=nodes,
#    mode="markers+text",
#    hoverinfo="text",
#    marker=dict(
#        color=[],
#        size=10,
#        line=dict(width=2)
#    )
#)
"""
fig = go.Figure(data=edge_trace + [node_trace],
                layout=go.Layout(
                    showlegend=False,
                    hovermode="closest",
                    margin=dict(b=20, l=5, r=5, t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    height=600,
                    width=800
                ))

"""