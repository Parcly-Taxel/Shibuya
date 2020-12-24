"""
The graphs here represent skeletons of polyhedra. If the graph is cubic then it likely has
a unit-distance embedding; if there are multiple higher-degree vertices it likely doesn't.
"""
from mpmath import sqrt, root, unitroots, polyroots
from shibuya.generators import cu, ring_edges, star_radius

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

def truncated_tetrahedron():
    """Return a unit-distance embedding of the truncated tetrahedron graph."""
    u4 = unitroots(4)
    z0 = sqrt(polyroots([1, -4+4j, -4j, 4+4j, -1])[1]) / 2
    r0 = [z0*u for u in u4]
    z1 = cu((-z0*1j).conjugate(), z0)
    r1 = [z1*u for u in u4]
    z2 = cu(-z0.conjugate(), z0)
    r2 = [z2*u for u in u4]
    vertices = r0 + r1 + r2
    edges = ring_edges(4, ((0, 0, 2), (0, 1, 0), (0, 2, 0), (1, 2, 0), (1, 2, -1)))
    return (vertices, edges)

def truncated_cube():
    """Return a unit-distance embedding of the truncated cube graph."""
    s1 = root(1, 12, 1)
    s2 = s1 + 1
    s3 = s1 + root(1, 6, 1)
    s4 = cu(s1, 0, 1, star_radius(3))
    u6 = unitroots(6)
    r1 = [s1*u for u in u6]
    r2 = [s2*u for u in u6]
    r3 = [s3*u for u in u6]
    r4 = [s4*u for u in u6]
    vertices = r1 + r2 + r3 + r4
    edges = ring_edges(6, ((0, 1, 0), (0, 2, 0), (1, 2, 0), (2, 1, 1), (3, 3, 2), (0, 3, 0)))
    return (vertices, edges)

def bidiakis_cube():
    """Return a unit-distance embedding of the bidiakis cube."""
    w6 = root(1, 6, 1)
    w8 = root(1, 8, 1)
    house = (0, 1j, w6, w6+1j, 1, 1+1j)
    house_edges = [(0, 1), (2, 3), (4, 5), (0, 2), (2, 4), (1, 3), (3, 5)]
    left_house = [v*w8 for v in house]
    right_house = [(v-1)*w8.conjugate() + 1 for v in house]
    right_house_edges = [(a + 6, b + 6) for (a, b) in house_edges]
    vertices = left_house + right_house
    edges = house_edges + right_house_edges + [(1, 6), (4, 11), (5, 7), (0, 10)]
    return (vertices, edges)
