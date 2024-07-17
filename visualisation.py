import rdflib
import plotly.graph_objects as go
import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go

# Parse the N-Triples file
g = rdflib.Graph()
g.parse("data.nt", format="nt")

# Extract the nodes (subjects and objects) and edges (predicates)
nodes = set()
edges = []
for s, p, o in g:
    nodes.add(str(s))
    nodes.add(str(o))
    edges.append((str(s), str(p), str(o)))

# Create the network graph
edge_trace = go.Scatter(
    x=[],
    y=[],
    mode="lines",
    line=dict(color="rgb(210,210,210)", width=1),
    hoverinfo="none"
)

for edge in edges:
    x0, y0 = nodes.index(edge[0]), nodes.index(edge[2])
    x1, y1 = nodes.index(edge[2]), nodes.index(edge[0])
    edge_trace["x"] += [x0, x1, None]
    edge_trace["y"] += [y0, y1, None]

node_trace = go.Scatter(
    x=[],
    y=[],
    text=[],
    mode="markers+text",
    hoverinfo="text",
    marker=dict(
        color=[],
        size=10,
        line=dict(width=2)
    )
)

for node in nodes:
    x, y = nodes.index(node), nodes.index(node)
    node_trace["x"].append(x)
    node_trace["y"].append(y)
    node_trace["text"].append(node)

fig = go.Figure(data=[edge_trace, node_trace],
                layout=go.Layout(
                    showlegend=False,
                    hovermode="closest",
                    margin=dict(b=20, l=5, r=5, t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    height=600,
                    width=800
                ))

#fig.show()


