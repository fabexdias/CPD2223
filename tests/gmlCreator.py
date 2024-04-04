import igraph as ig

# Read input file
with open('gen19-23.in', 'r') as f:
    nb_nodes, nb_edges = map(int, f.readline().split())
    edges = [tuple(map(int, line.split()[:2])) for line in f]

# Create graph
g = ig.Graph()
g.add_vertices(nb_nodes)
g.add_edges(edges)

# Write GML file
ig.write(g, 'gen19-23.gml', format='gml')
