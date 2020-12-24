"""
Other unit-distance graphs...
"""
from mpmath import mpc, sqrt, root, unitroots
from shibuya.generators import cu, star_radius, ring_edges

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
