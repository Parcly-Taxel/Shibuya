"""
These functions generate simple collections of vertices or edges, or make new ones from old.
"""
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

def cartesian_product(*graphs):
    """Given a list of graphs, construct their Cartesian product.
    The metrics of each graph are preserved, so if they are all unit-distance
    the product is also unit-distance."""
    pass

def cu(z1, z2):
    """Constructs the point at a distance of 1 (hence Construct Unit) from z1 and z2,
    left of the line from z1 to z2."""
    m = (z1 + z2) / 2
    d, theta = polar(z2 - z1)
    return m + sqrt(1 - d * d / 4) * expj(theta + pi / 2)

def cc(z1, z2, r1, r2):
    """Constructs the point at a distance of r1 from z1 and r2 from z2, left of the
    line from z1 to z2 (hence Construct Cosine (rule))."""
    d, theta0 = polar(z2 - z1)
    theta = acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1))
    return z1 + r1 * expj(theta0 + theta)

def subtend_r(p, q):
    """Calculates the radius of the circle if a unit-length chord on it subtends
    p/q of a full rotation. This is used to precisely place points on a star of edges."""
    return sqrt(0.5 / (1 - cospi(2 * mpf(p) / mpf(q))))

def ring_edges(N, triples):
    """Suppose the vertices are grouped into some number of rings, each ring having N vertices.
    Those vertices can then be grouped into N congruent units.
    This function returns an index list for the edges constructed from triples,
    where a triple (a, b, k) means "a vertex in ring a is connected to the vertex in ring b, k units forward."."""
    L = []
    for (a, b, k) in triples:
        L.extend([(i + a * N, (i + k) % N + b * N) for i in range(N)])
    return L
