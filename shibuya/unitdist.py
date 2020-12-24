from math import gcd
from mpmath import root, unitroots
from shibuya.generators import *

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

# Most of the following embeddings were taken from MathWorld

def pappus():
    """Return a unit-distance embedding of the Pappus graph (F18A)."""
    u6 = unitroots(6)
    r0 = [u*0.5j for u in u6]
    z1 = cu(r0[2], r0[0])
    r1 = [z1*u for u in u6]
    z2 = cu(0, z1)
    r2 = [z2*u for u in u6]
    vertices = r0 + r1 + r2
    edges = ring_edges(6, ((0, 0, 3), (0, 1, 0), (0, 1, -2), (2, 2, 1), (1, 2, 0)))
    return (vertices, edges)

def coxeter():
    """Return a unit-distance embedding of the Coxeter graph (F28A)."""
    u7 = unitroots(7)
    s1 = star_radius(7)
    s2 = star_radius(7, 2)
    s3 = star_radius(7, 3)
    r0 = [-s2*u for u in u7]
    r1 = [s3*u for u in u7]
    z2 = cu(r0[0], r1[3])
    r2 = [z2*u for u in u7]
    z3 = cu(0, z2, s1, 1)
    r3 = [z3*u for u in u7]
    vertices = r0 + r1 + r2 + r3
    edges = ring_edges(7, ((0, 0, 2), (1, 1, 3), (3, 3, 1), (0, 2, 0), (1, 2, -3), (2, 3, 0)))
    return (vertices, edges)

def tutte8_vertices(x):
    u5 = unitroots(5)
    r0 = [u*0.25j for u in u5]
    r1 = [-u*x*1j for u in u5]
    z2 = cu(r0[2], r0[4])
    r2 = [z2*u for u in u5]
    z3 = cu(r0[3], r1[1])
    r3 = [z3*u for u in u5]
    z4 = cu(r1[4], r1[3])
    r4 = [z4*u for u in u5]
    z5 = cu(r3[0], r2[0])
    r5 = [z5*u for u in u5]
    return (r0 + r1 + r2 + r3 + r4 + r5, abs(z4-z5)-1)

def tutte8():
    """Return a unit-distance embedding of the Tutte 8-cage (F30A; MathWorld calls this
    the Levi graph)."""
    x0 = findroot(lambda x: tutte8_vertices(x)[1], [0.33, 0.34])
    vertices = tutte8_vertices(x0)[0]
    edges = ring_edges(5, ((0, 2, 1), (0, 2, 3), (0, 3, 2), (1, 3, 4), (1, 4, 1), (1, 4, 2),
                           (2, 5, 0), (3, 5, 0), (4, 5, 0)))
    return (vertices, edges)
