"""
Braced polygons, triangle-free and otherwise.
All graphs are Laman unless otherwise stated; if such a graph has v vertices
it has 2v-3 edges.
"""
from mpmath import *
from shibuya.generators import cu, star_radius, ring_edges, all_unit_distances, delete_vertices
from shibuya.rigidity import jacobian

def khodulyov_square():
    """Return Andrei Khodulyov's 11-vertex rigid square."""
    s1 = star_radius(4)
    r1 = [s1*u for u in unitroots(4)]
    r1.extend([r1[1]-1, r1[2]-1, r1[3]-1])
    r1.extend([cu(r1[1], r1[6]), cu(r1[6], r1[1])])
    r1.extend([cu(r1[3], r1[4]), cu(r1[4], r1[3])])
    return all_unit_distances(r1)

def khodulyov_pentagon():
    """Return Andrei Khodulyov's 17-vertex rigid regular pentagon."""
    A = 0
    B = 1
    p0 = expj(-2*pi/5)
    p1 = 1+expj(-3*pi/5)
    p2 = cu(p1, p0)
    p3 = cu(p2, A)
    p4 = cu(B, p2)
    p5 = cu(p3, p1)
    p6 = cu(p0, p4)
    p7 = cu(p1, A)
    p8 = cu(p7, p5)
    p9 = cu(p6, B)
    p10 = cu(p4, p8)
    p11 = cu(p8, p4)
    p12 = cu(p9, p8)
    p13 = cu(p8, p9)
    p14 = cu(p8, B)
    vertices = (A, B, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14)
    return all_unit_distances(vertices)

def rigid_heptagon(suppress=None):
    """Return a unit-distance graph containing a regular heptagon that was conjectured by
    Ed Pegg to be rigid in https://math.stackexchange.com/q/3954719/357390.
    I proved it to be rigid, even when two adjacent vertices are deleted;
    use (2,3) or (2,2) or (3,3) as the parameter to see these reduced rigid graphs
    which all have 19 vertices."""
    u7 = unitroots(7)
    s1 = star_radius(7)
    r1 = [s1*u for u in u7]
    z2 = cu(0, s1, star_radius(7, 2), 1)
    r2 = [z2*u for u in u7]
    z3 = cu(s1, 0, 1, star_radius(7, 3))
    r3 = [z3*u for u in u7]
    vertices = r1 + r2 + r3
    edges = ring_edges(7, ((0, 0, 1), (1, 1, 2), (0, 1, 0), (2, 2, 3), (0, 2, 0), (1, 2, 0)))
    G = (vertices, edges)
    if suppress == (2,3):
        return delete_vertices(G, (7, 14))
    if suppress == (2,2):
        return delete_vertices(G, (8, 13))
    if suppress == (3,3):
        return delete_vertices(G, (14, 17))
    return G

def rigid_14gon():
    """Return Somsky's braced 14-gon (47 vertices), based on the rigid heptagon."""
    core = rigid_heptagon((2, 3))[0]
    core = [v - core[0] for v in core]
    z3, z2, z1 = core[6:3:-1]
    r3 = [z3*root(1, 14, n+1) for n in range(7)]
    r2 = [z2*root(1, 14, n+1) for n in range(9)]
    r1 = [z1*root(1, 14, n+1) for n in range(12)]
    vertices = core + r3 + r2 + r1
    return all_unit_distances(vertices)

def cloud9_vertices(t, u):
    p0 = 0
    p1 = 1
    p2 = mpc(t, sqrt(1-t**2))
    p3 = mpc(sqrt(1-u**2), u)
    p4 = p2 + 1
    p5 = p3 + 1
    p6 = cu(p3, p4)
    p7 = cu(p4, p3)
    p8 = cu(p2, p5)
    p9 = cu(p8, p6)
    p10 = cu(p7, p8)
    p11 = cu(p10, p6)
    vertices = (p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11)
    return vertices, abs(p0 - p11) - 1, abs(p1 - p9) - 1

def cloud9():
    """Minimal triangle-free braced square, 12 vertices, based on the hodfish.
    Looks like Cloud9's logo if you squint."""
    f = lambda *x: cloud9_vertices(*x)[1:]
    x0 = (0, 0.75)
    print(f(*x0)) # zero within floating-point error
    print(det(jacobian(f, x0))) # non-zero
    vertices = cloud9_vertices(*x0)[0]
    return all_unit_distances(vertices)

def tfrhexagon_vertices(v, mode=0):
    p0 = 0
    p1 = -1
    p2 = root(1, 6, 1)
    p3 = root(1, 6, 5)
    p4 = p1 + p2 # v1
    p5 = p3 + p1 # v2
    p6 = p2 + p3
    p7 = expj(v) # v3
    p8 = p4 + p7
    p9 = p5 + p7
    p10 = cu(p8, p5) if mode == 1 else cu(p5, p8)
    p11 = p10 - p5 # v4
    p12 = p7 + p11
    p13 = p4 + p7 + p11
    p14 = p4 + p11
    p15 = p5 + p7 + p11
    vertices = (p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15)
    return vertices, abs(p15 - p6) - 1

def tfrhexagon(mode=0):
    """Triangle-free braced hexagon in 16 vertices. mode (0 or 1) selects between
    two algebraically related versions (corresponding coordinates have the same
    minimal polynomial)."""
    f = lambda v: tfrhexagon_vertices(v, mode)[1]
    v0 = findroot(f, -0.6)
    print(diff(f, v0)) # non-zero
    vertices = tfrhexagon_vertices(v0, mode)[0]
    edges = list(all_unit_distances(vertices)[1])
    edges.remove((0, 4))
    edges.remove((0, 5))
    edges.remove((0, 6))
    return (vertices, edges)
