"""
These functions generate simple collections of vertices or edges, or make new ones from old.
"""
from functools import reduce
from mpmath import *

def disjoint_union(*graphs):
    """Given a list of graphs, construct their disjoint union."""
    res_vertices = []
    res_edges = []
    for (vertices, edges) in graphs:
        l = len(res_vertices)
        res_edges.extend((a+l, b+l) for (a, b) in edges)
        res_vertices.extend(vertices)
    return (res_vertices, res_edges)

def cartesian_product2(g1, g2):
    """Given two graphs, construct their Cartesian product.
    The metrics of each graph are preserved, so if they are all unit-distance
    the product is also unit-distance.

    g1 is held steady and the first vertex of g2 is put on each of g1's vertices
    to generate all the product's vertices."""
    vertices1, edges1 = g1
    vertices2, edges2 = g2
    anchor = vertices2[0]
    n = len(vertices1)
    vertices = [v1 + v2 - anchor for v1 in vertices1 for v2 in vertices2]
    edges = []
    for (a, b) in edges1:
        m = len(vertices2)
        for i in range(m):
            edges.append((a*m+i, b*m+i))
    for (c, d) in edges2:
        n = len(vertices1)
        for j in range(n):
            edges.append((j*m+c, j*m+d))
    return (vertices, edges)

def cartesian_product(*graphs):
    """Given a list of graphs, construct their Cartesian product using the
    conventions in cartesian_product2()."""
    return reduce(cartesian_product2, graphs)

def cu(z1, z2, r1=1, r2=1):
    """Constructs the point at a distance of r1 from z1 and r2 from z2, left of the
    line from z1 to z2. If r1 or r2 is omitted the default is 1."""
    d, theta0 = polar(z2 - z1)
    if not (abs(r1-r2) <= d <= r1+r2):
        raise ValueError("point cannot be constructed")
    theta = acos((r1*r1 - r2*r2 + d*d) / (2*d*r1))
    return z1 + r1 * expj(theta0 + theta)

def star_radius(p, q=1):
    """Calculates the radius of the regular star polygon {p/q} with
    unit-length edges. If q=1 the star becomes a polygon."""
    return 0.5 / sinpi(fraction(q,p))

def ring_edges(n, triples):
    """Let there be n*r vertices grouped into r ordered rings of n vertices each.
    Each triple (a, b, k) in triples corresponds to the statement "for 0 <= i <= n, vertex i
    in ring a is connected to vertex i+k in ring b, all indices taken modulo n".
    Return all the edges implied by these statements."""
    L = []
    for (a, b, k) in triples:
        limit = k if a == b and 2*k == n else n
        L.extend([(i+a*n, (i+k)%n + b*n) for i in range(limit)])
    return L

def all_unit_distances(vertices, tol=1e-12):
    """Returns the graph formed by inserting all edges of length 1 between the vertices."""
    n = len(vertices)
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            if almosteq(abs(vertices[i]-vertices[j]), 1, tol):
                edges.append((i, j))
    return (vertices, edges)

def delete_vertices(graph, dverts):
    """Delete the vertices indexed by dverts from graph; return the resulting graph."""
    vertices, edges = graph
    remverts = list(filter(lambda x: x not in dverts, range(len(vertices))))
    vmap = {v: n for (n, v) in enumerate(remverts)}
    nvertices = list(map(vertices.__getitem__, remverts))
    nedges = []
    for (a, b) in edges:
        if a in dverts or b in dverts:
            continue
        nedges.append((vmap[a], vmap[b]))
    return (nvertices, tuple(nedges))
