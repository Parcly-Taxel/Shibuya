"""
Cubic symmetric graphs with more than 102 vertices.
"""
from mpmath import *
from shibuya.generators import star_radius
from shibuya.generators import fixparams_unitdist, symmetrise, remove_edges

@fixparams_unitdist(0.37, 1.47)
def f104a(a, b):
    p1 = mpc(a, 0.5)
    p2 = mpc(b, 0.5)
    d1 = abs(p1 - p2*root(1, 26, 4)) - 1
    d2 = abs(p2 - p1/root(1, 26, 1)) - 1
    return (symmetrise((p1, p2), "D26"), (d1, d2))

@remove_edges(lambda e: {e[0]//9, e[1]//9} < {0, 2, 4} and e[0]//9 != e[1]//9)
@fixparams_unitdist(0.18, 3, 3)
def f108a(a, b, c):
    z1 = -star_radius(9)
    z2 = z1 - 1
    z3 = -star_radius(9, 2)
    z4 = z3 + 1
    z5 = star_radius(9, 4)
    z6 = z5 + 1
    z7 = z2 + expj(a)
    z8 = z4 + expj(b)
    z9 = z6 + expj(c)
    d1 = abs(z7 - z8*root(1, 9, 1)) - 1
    d2 = abs(z8 - conj(z9)*root(1, 9, -2)) - 1
    d3 = abs(z9 - conj(z7)*root(1, 9, 5)) - 1
    vertices = symmetrise((z1, z2, z3, z4, z5, z6), "C9") + symmetrise((z7, z8, z9), "D9")
    return (vertices, (d1, d2, d3))
