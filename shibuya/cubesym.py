"""
Cubic symmetric graphs. Most of the embeddings realised here were taken from MathWorld.
"""
from mpmath import sqrt, findroot, unitroots
from shibuya.generators import cu, star_radius, ring_edges

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

def dyck():
    """Return a unit-distance embedding of the Dyck graph (F32A)."""
    r0 = unitroots(8)
    r1 = [sqrt(2)*u for u in r0]
    z2 = cu(r0[1], 0, 1, star_radius(8))
    r2 = [z2*u for u in r0]
    z3 = cu(0, r1[0], star_radius(8, 3), 1)
    r3 = [z3*u for u in r0]
    vertices = r0 + r1 + r2 + r3
    edges = ring_edges(8, ((0, 1, 1), (0, 1, -1), (0, 2, -1), (2, 2, 1), (1, 3, 0), (3, 3, 3)))
    return (vertices, edges)

def biggssmith():
    """Return a unit-distance embedding of the Biggs–Smith graph (F102A)."""
    s1 = star_radius(17)
    s2 = star_radius(17, 2)
    s4 = star_radius(17, 4)
    s8 = star_radius(17, 8)
    u17 = unitroots(17)
    r1 = [s1*u*1j for u in u17]
    r4 = [s4*u*1j for u in u17]
    r8 = [s8*u*-1j for u in u17]
    sh1 = cu(r1[0], r4[0])
    rh1 = [sh1*u for u in u17]
    sh2 = cu(sh1, r8[7])
    rh2 = [sh2*u for u in u17]
    s2 = cu(sh2, 0, 1, s2)
    r2 = [s2*u for u in u17]
    vertices = r1 + r4 + rh1 + r8 + rh2 + r2
    edges = ring_edges(17, ((0, 0, 1), (1, 1, 4), (3, 3, 8), (5, 5, 2),
                            (0, 2, 0), (1, 2, 0), (2, 4, 0), (4, 5, 0), (4, 3, 7)))
    return (vertices, edges)
