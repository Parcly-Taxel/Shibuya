"""
The I-graphs are always unit-distance. The functions here generate such embeddings
and include a convenience function for generalised Petersen graphs.
"""
from math import gcd
from mpmath import sqrt, root, unitroots, almosteq
from shibuya.generators import disjoint_union, cartesian_product2, cu, star_radius, ring_edges

def igraph(n, j, k):
    """Return a unit-distance embedding of the I-graph (n,j,k), either directly
    or through an isomorph. Based on Žitnik, Horvat and Pisanski (2010),
    All generalized Petersen graphs are unit-distance graphs,
    http://preprinti.imfm.si/PDF/01109.pdf"""
    if (d := gcd(n,gcd(j,k))) > 1:
        verts, edges = igraph(n//d, j//d, k//d)
        vgroups = [[v * root(1, n, a) for v in verts] for a in range(d)]
        return disjoint_union(*((vs, edges) for vs in vgroups))
    j = min(j % n, -j % n)
    k = min(k % n, -k % n)
    if j == k:
        r = star_radius(n)
        cycle = [r*u + (0.5j if n % 4 == 2 else 0.5) for u in unitroots(n)]
        edge = (0, -1j if n % 4 == 2 else -1)
        cycle_edges = ring_edges(n, [(0, 0, 1)])
        return cartesian_product2((cycle, cycle_edges), (edge, [(0, 1)]))
    j, k = sorted((j, k))
    if (n, j, k) == (12, 1, 5): # Nauru graph case
        R = (6 + 7*sqrt(3) + sqrt(15)) / 12
        r = (6 + 7*sqrt(3) - sqrt(15)) / 12
        L = R-1
        l = r-1
        u12 = unitroots(12)
        ring_R = [u*R for u in u12[::2]]
        ring_r = [u*r for u in u12[1::2]]
        ring_L = [u*L for u in u12[::2]]
        ring_l = [u*l for u in u12[1::2]]
        vertices = ring_R + ring_r + ring_L + ring_l
        edges = ring_edges(6, ((0, 1, 0), (0, 1, -1), (0, 2, 0), (1, 3, 0), (2, 3, 2), (2, 3, -3)))
        return (vertices, edges)
    # Find an admissible isomorph
    for a in range(1, n):
        if gcd(n, a) > 1:
            continue
        aj = a*j % n
        ak = a*k % n
        r1 = star_radius(n, aj)
        r2 = star_radius(n, ak)
        dist = abs(r1 - r2)
        if dist < 1 or almosteq(dist, 1):
            z2 = cu(0, r1, r2, 1)
            vertices = [r1*u for u in unitroots(n)] + [z2*u for u in unitroots(n)]
            edges = ring_edges(n, ((0, 0, aj), (1, 1, ak), (0, 1, 0)))
            return (vertices, edges)
    raise ValueError("can't happen")

named_gp = {"cube": (4, 1),
            "petersen": (5, 2),
            "durer": (6, 2),
            "mobiuskantor": (8, 3),
            "dodecahedron": (10, 2),
            "desargues": (10, 3),
            "nauru": (12, 5)}

def genpetersen(n=5, k=2):
    """Return a unit-distance embedding of the generalised Petersen graph GP(n, k),
    equivalent to I(n, 1, k). n can be a string naming a specific graph of this family."""
    if type(n) == str:
        n, k = named_gp[n]
    return igraph(n, 1, k)
