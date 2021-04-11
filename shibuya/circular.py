"""
Functions to draw graphs with a distinctly circular layout.
They are not always unit-distance or integral, obviously.
"""
from mpmath import root, unitroots
from shibuya.generators import star_radius, ring_edges, lcf_edges, cartesian_product

def circulant(n, taps):
    """Return the circulant graph on n vertices with offsets given by taps.
    The vertices are placed on an n-gon of unit side length."""
    r = star_radius(n)
    vertices = [r*u for u in unitroots(n)]
    edges = ring_edges(n, [(0, 0, k) for k in taps])
    return (vertices, edges)

def cycle(n):
    """Return the cycle graph on n vertices."""
    return circulant(n, [1])

def complete(n):
    """Return the complete graph on n vertices."""
    return circulant(n, range(1, n//2+1))

def lcf_graph(n, pattern):
    """Return the graph with the given LCF notation; the number of repeats is implied.
    Vertices are arranged in a circle."""
    r = star_radius(n)
    vertices = [r*u for u in unitroots(n)]
    return (vertices, lcf_edges(n, pattern))

def mobiusladder(n=3):
    """Return the Möbius ladder on 2n vertices. The n = 3 case corresponds
    to the utility graph; its embedding returned by this function is integral
    with edge lengths 1 and 2, as well as rigid, but not first-order rigid."""
    return circulant(2*n, (1, n))

def hamming(d, q):
    """Return the Hamming graph H(d, q), the Cartesian product of d copies of K_q.
    H(d, 2) is the hypercube graph Q_d. H(d, 2) and H(d, 3) are unit-distance
    (and are rendered in such a fashion by this function), but have d-1 degrees
    of freedom."""
    verts, edges = complete(q)
    vgroups = [[v * root(1, d*q, a) for v in verts] for a in range(d)]
    return cartesian_product(*((vs, edges) for vs in vgroups))
