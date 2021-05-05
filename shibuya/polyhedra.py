"""
Skeletons of polyhedra. If the graph is cubic then it likely has a unit-distance embedding;
if there are multiple higher-degree vertices it likely doesn't.
"""
from mpmath import *
from shibuya.generators import cu, ring_edges, star_radius
from shibuya.generators import fixparams_unitdist, symmetrise, remove_edges

def tetrahedron():
    """Return the integral embedding of the tetrahedral graph with shortest
    longest edge. This graph, octahedron() and icosahedron() are from
    Harborth and Moller, Minimum Integral Drawings of the Platonic Graphs,
    Mathematics Magazine vol. 67 no. 5 (December 1994), 355-358."""
    v1 = 0
    v2 = 4
    v3 = cu(v1, v2, 2, 4)
    v4 = cu(v1, v2, 4, 2)
    vertices = (v1, v2, v3, v4)
    edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    return (vertices, edges)

def octahedron():
    """Return the best integral embedding of the octahedral graph."""
    v1 = (13 * root(1, 6, 1) - 2) / 3
    v2 = -v1.conjugate()
    omega = root(1, 3, 1)
    v3 = v1 * omega
    v4 = v2 * omega
    v5 = v1 / omega
    v6 = v2 / omega
    vertices = (v1, v2, v3, v4, v5, v6)
    edges = ((0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3),
             (1, 5), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5))
    return (vertices, edges)

def icosahedron():
    """Return the best integral embedding of the icosahedral graph."""
    v0 = 0
    v1 = 8
    v2 = cu(v0, v1, 4, 8)
    v3 = v2 + 8
    v4 = cu(v0, v1, 7, 6)
    v5 = v4 + 8
    v7 = cu(v2, v3, 4, 6)
    v6 = v7 - 8
    v9 = cu(v4, v5, 4, 6)
    v8 = v9 - 8
    v10 = cu(v6, v7, 7, 6)
    v11 = v10 + 8
    vertices = (v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11)
    edges = ((0, 1), (0, 2), (0, 4), (0, 6), (0, 8), (1, 2),
             (1, 3), (1, 4), (1, 5), (2, 3), (2, 6), (2, 7),
             (3, 5), (3, 7), (3, 11), (4, 5), (4, 8), (4, 9),
             (5, 9), (5, 11), (6, 7), (6, 8), (6, 10), (7, 10),
             (7, 11), (8, 9), (8, 10), (9, 10), (9, 11), (10, 11))
    return (vertices, edges)

@fixparams_unitdist()
def truntet():
    """Return a unit-distance embedding of the truncated tetrahedron graph."""
    b = sqrt(sqrt(3)-sqrt(4*sqrt(3)-5)) / 2
    z0 = mpc(sqrt(1-b**2), b) / 2
    z1 = cu(conj(z0)*1j, z0)
    z2 = cu(-conj(z0), z0)
    return symmetrise((z0, z1, z2), "C4")

@remove_edges(lambda e: e[0]//6 == e[1]//6 == 2)
@fixparams_unitdist()
def truncube():
    """Return a unit-distance embedding of the truncated cube graph."""
    z0 = 1j + root(1,6,1)
    z1 = 1j + root(1,3,1)
    z3 = cu(1j, 0, 1, star_radius(3))
    return symmetrise((z0, z1, 1j, z3), "C6")

@remove_edges(lambda e: e[0]//6 == e[1]//6 == 0)
@fixparams_unitdist()
def trunoct():
    """Return a unit-distance embedding of the truncated octahedron graph."""
    z0 = star_radius(3)*1j
    z1 = 2*z0
    a = (sqrt(33)-3) / 12
    z2 = z0 + mpc(a, sqrt(1-a**2))
    z3 = z2 + 1
    return symmetrise((z0, z1, z2, z3), "C6")

@fixparams_unitdist(-2.28)
def trundodec(a):
    """Return a unit-distance embedding of the truncated dodecahedron graph."""
    p0 = star_radius(20)*root(1,40,1)
    p1 = star_radius(20)*root(1,40,3)
    p2 = cu(p0, p1)
    p3 = p2 + expj(a)
    p4 = p3 - 1
    p5 = cu(p3, p4)
    return (symmetrise((p0, p1, p2, p3, p4, p5), "D5"), [abs(p4 - root(1,5,1)*p5) - 1])

@fixparams_unitdist(3.22)
def trunicos(b):
    """Return a unit-distance embedding of the truncated icosahedron graph."""
    p0 = star_radius(5)*root(1,20,1)
    p1 = p0 + root(1,20,1)
    p2 = mpc(b, 0.5)
    p3 = cu(p2, p1)
    p4 = cu(p3, p1*root(1,5,-1))
    p5 = cu(p4, p2*root(1,5,-1))
    return (symmetrise((p0, p1, p2, p3, p4, p5), "D5"),
            [abs(p5 - root(1,5,-1)*conj(p5)) - 1])

@fixparams_unitdist()
def bidiakis():
    """Return a unit-distance embedding of the bidiakis cube."""
    a = (-sqrt(7)+3j) / 4
    return symmetrise((0.5, a+0.5, a-0.5), "C4")
