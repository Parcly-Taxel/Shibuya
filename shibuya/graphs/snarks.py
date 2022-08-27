"""
Unit-distance embeddings of snarks (3-connected cubic graphs with chromatic index 4),
including the Tietze graph.
"""
from mpmath import *
from shibuya.generators import (cu, star_radius, ring_edges,
        all_unit_distances, fixparams_unitdist, symmetrise)

def tietze():
    """Return a unit-distance embedding of Tietze's graph."""
    t = (13*sqrt(3) + sqrt(13 * (40*sqrt(3)-9))) / 52
    b = 1 / sqrt(3)
    c = b + root(1, 3, 1)
    c2 = c * root(1, 3, 1)
    d = c - mpc(t, sqrt(1-t*t))
    e = cu(d, c2)
    vertices = [v * u for u in unitroots(3) for v in (b, c, d, e)]
    edges = [(0, 1), (1, 2), (2, 3), (0, 4), (3, 5), (2, 7),
             (4, 5), (5, 6), (6, 7), (4, 8), (7, 9), (6, 11),
             (8, 9), (9, 10), (10, 11), (8, 0), (11, 1), (10, 3)]
    return (vertices, edges)

def blanusa1():
    """Draws a unit-distance embedding of the first Blanuša snark."""
    p0 = rect(0.5, atan(1/sqrt(2)))
    p1 = -conj(p0)
    p2 = -p0
    p3 = conj(p0)
    p4 = cu(p1, p0)
    p5 = cu(p2, p1)
    p6 = cu(p3, p2)
    p7 = cu(p0, p3)
    s1 = [p0, p1, p2, p3, p4, p5, p6, p7]
    s2 = [p-1j for p in s1]
    A = cu(s1[2], s1[0])
    B = cu(s1[1], s1[3])
    vertices = s1 + s2 + [A, B]
    edges = set(all_unit_distances(vertices)[1])
    edges -= {(0, 8), (1, 9), (2, 10), (3, 11), (0, 2), (1, 3)}
    return (vertices, list(sorted(edges)))

def blanusa2_vertices(t, u):
    p0 = rect(0.5, t)
    p1 = -conj(p0)
    p2 = -p0
    p3 = conj(p0)
    p4 = cu(p1, p0)
    p5 = cu(p2, p1)
    p6 = cu(p3, p2)
    p7 = cu(p0, p3)
    s1 = [p0, p1, p2, p3, p4, p5, p6, p7]
    s2 = [p+expj(u) for p in s1]
    A = cu(s2[4], s1[4])
    B = cu(s1[7], s2[7])
    vertices = s1 + s2 + [A, B]
    return vertices, abs(A - B) - 1

def blanusa2(t=pi/3):
    """Draws a unit-distance embedding of the second Blanuša snark.
    t (0 <= t <= pi/2) controls the proportions of the two stars inside."""
    f = lambda u: blanusa2_vertices(t, u)[1]
    u0 = findroot(f, 0.1)
    vertices = blanusa2_vertices(t, u0)[0]
    edges = set(all_unit_distances(vertices)[1])
    edges -= {(0, 8), (1, 9), (2, 10), (3, 11), (4, 12), (7, 15)}
    return (vertices, list(sorted(edges)))

def flowersnark(n=5):
    """Return a unit-distance embedding of the flower snark J_n,
    where n is an odd number at least 5."""
    un = unitroots(n)
    s0 = 2*star_radius(n, 2)
    r0 = [s0*u for u in un]
    s1 = r0[1].real
    r1 = [s1*u for u in un]
    z2 = cu(s1, s0)
    r2 = [z2*u for u in un]
    z3 = cu(0, z2, star_radius(n), 1)
    r3 = [z3*u for u in un]
    vertices = r0 + r1 + r2 + r3
    edges = ring_edges(n, ((1, 0, 1), (1, 0, -1), (0, 2, 0), (1, 2, 0), (2, 3, 0), (3, 3, 1)))
    return (vertices, edges)

@fixparams_unitdist(-1.2, -0.375, 2.35, 0.5)
def szekeres(a, b, c, d):
    """Return a unit-distance embedding of the Szekeres snark."""
    u = root(1,5,1)
    p1 = a+1j
    p2 = b+0.5j
    p3 = c+0.5j
    p4 = a-1 + expj(d)
    cons = (abs(p4 - u*u*p3) - 1, abs(p4 - u*p2) - 1,
            abs(p1 - u*u*p3) - 1, abs(p1 - conj(u*p2)) - 1)
    return (symmetrise((a, a-1), "C5") + symmetrise((p1, p2, p3, p4), "D5"), cons)

@fixparams_unitdist(0.7)
def watkins(a):
    """Return a unit-distance embedding of the Watkins snark."""
    p0 = 0.49-0.32j
    p1 = (p0-1) * root(1,5,2)
    s0 = p0 + expj(0.06)
    s1 = s0 + expj(0.02)
    s2 = cu(p1, s1)
    s3 = s0 + expj(a)
    p2 = cu(p1, s3)
    p3 = cu(p2, s1)
    s4 = cu(s3, s2)
    p4 = cu(s4, p0)
    return (symmetrise((p0, p1, p2, p3, p4, s0, s1, s2, s3, s4), "C5"), [abs(p3 - root(1,5,1)*p4) - 1])
