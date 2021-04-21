"""
Cubic symmetric graphs. Most of the embeddings realised here were taken from MathWorld.
"""
from mpmath import *
from functools import reduce
from shibuya.generators import cu, star_radius, ring_edges, lcf_edges
from shibuya.generators import all_unit_distances, circumcentre
from shibuya.generators import fixparams_unitdist, symmetrise

# F4A = tetrahedron() or complete(4) (not unit-distance)
# F6A = circulant(6, (1, 3)) or mobiusladder(3) (not unit-distance)
# F8A = genpetersen("cube")
# F10A = genpetersen("petersen")

def heawood():
    """Return the symmetric unit-distance embedding of the Heawood graph (F14A)
    tucked away in Mathematica's GraphData."""
    P = [10485760, 78643200, 263192576, 543686656, 812777472, 942080000, 843317248, 552468480, 208879616, -31170560, -99213312, -76779520, -32795648, 7878144, 17269760, 16256512, 11392032, 4836080, 3014064, 361320, 69498, -165789]
    v0 = findroot(lambda v: polyval(P, v), 0.275)
    p0 = mpc(0.5, v0)
    p1 = mpc(sqrt(1-(v0+0.5)**2)-0.5, -0.5)
    p2 = cu(p0, -p0)
    p3 = cu(p1, -p1)
    p4 = cu(p2, p3)
    vertices = [mpc(s*re(v), im(v)) for s in (1, -1) for v in (p0, -p0, p1, -p1, p2, p3, p4)]
    return all_unit_distances(vertices)

# F16A = genpetersen("mobiuskantor")

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

# F20A = genpetersen("dodecahedron")
# F20B = genpetersen("desargues")
# F24A = genpetersen("nauru")

def f26a_vertices(t):
    A, B, C = unitroots(6)[4:1:-1]
    p2 = mpc(t, sqrt(1-t**2)) / 2
    p1 = p2 * root(1, 6, 1)
    p3 = p2 * root(1, 6, 5)
    p4 = cu(p1, B)
    p5 = cu(p2, C)
    p6 = cu(p3, -A)
    p7 = cu(p1, -p6)
    p8 = cu(p4, p2)
    p9 = cu(p5, p3)
    p10 = cu(-p8, p7)
    V = (A, B, C, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
    return ([v*s for v in V for s in (1, -1)], abs(p9 - p10) - 1)

def f26a():
    """Return a unit-distance embedding of the F26A graph."""
    t0 = findroot(lambda t: f26a_vertices(t)[1], 0.2)
    return all_unit_distances(f26a_vertices(t0)[0])

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

@fixparams_unitdist(-1.76, -0.76)
def tutte8(a, b):
    z1 = 1.69j
    z2 = z1 + expj(a)
    z3 = z2 + sign(-z2)
    z4 = cu(z1*root(1,5,2), z1*root(1,5,-2))
    d1 = abs(z3-z4) - 1
    z5 = z2 + expj(b)
    z6 = cu(z5*root(1,5,2), z5)
    d2 = abs(z6-z3*root(1,5,1)) - 1
    vertices = symmetrise((z1, z2, z3, z4, z5, z6), "C5")
    return (vertices, (d1, d2))

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

def dyck2_vertices(a):
    p1 = mpc(a, 0.5)
    p2 = mpc(3*a, 0.5)
    p3 = cu(p1, -1j*p1)
    p4 = cu(p2, -1j*p2)
    vertices = [u*p for u in (1, 1j, -1, -1j) for p in (p1, p2, p3, p4)]
    vertices.extend([conj(p) for p in vertices])
    return (vertices, abs(p3 - conj(p4)) - 1)

def dyck2():
    """Return a unit-distance embedding of the Dyck graph with D4 symmetry."""
    t0 = findroot(lambda t: dyck2_vertices(t)[1], 0.1)
    return all_unit_distances(dyck2_vertices(t0)[0])

def f38a_vertices(t):
    u6 = unitroots(6)
    vertices = []
    p0 = 1
    p1 = 2
    p2 = rect(0.5, t)
    p3 = cu(p2, p1)
    p4 = cu(p3, 2*u6[1])
    for (i, u) in enumerate(u6):
        p5 = cu(u6[4]*p2, p4) if i in (2, 5) else cu(p4, u6[4]*p2)
        vertices.extend(p*u for p in (p0, p1, p2, p3, p4, p5))
    p6 = circumcentre(vertices[5], vertices[17], vertices[29])
    vertices.extend((p6, -p6))
    return (vertices, abs(p6 - vertices[5]) - 1)

def f38a():
    """Return a unit-distance embedding of the F38A graph."""
    t0 = findroot(lambda t: f38a_vertices(t)[1], 0.29)
    return all_unit_distances(f38a_vertices(t0)[0])

def f40a(x=0.75):
    """Return a unit-distance embedding of F40A (bipartite double cover of F20A).
    x can be anything between (sqrt(5)-1)/2 and 1."""
    u10 = unitroots(10)
    z0 = star_radius(10)
    r0 = [z0*u for u in u10]
    z1 = cu(r0[1], 0, 1, x)
    r1 = [z1*u for u in u10]
    z2 = cu(r1[2], r1[-2])
    r2 = [z2*u for u in u10]
    z3 = cu(0, z2, z0, 1)
    r3 = [z3*u for u in u10]
    vertices = r0 + r1 + r2 + r3
    return all_unit_distances(vertices)

def f42a_vertices(a, b, c):
    u7 = unitroots(7)
    pa = mpc(a, 0.5)
    pb = mpc(b, 0.5)
    pc = mpc(c, 0.5)
    pac, pbc, pcc = (conj(p) for p in (pa, pb, pc))
    d1 = abs(pa - u7[1]*pbc)**2 - 1
    d2 = abs(pb - u7[2]*pcc)**2 - 1
    d3 = abs(pc - u7[4]*pac)**2 - 1
    vertices = [u*p for u in u7 for p in (pa, pb, pc, pac, pbc, pcc)]
    return (vertices, (d1, d2, d3))

def f42a(mode=0):
    """Return a unit-distance embedding of the F42A graph.
    mode (0 or 1) selects between two algebraically related forms."""
    x0 = (0.27, 1.36, 0.52) if mode == 0 else (1.24, 0.18, -0.53)
    t0 = findroot(lambda *t: f42a_vertices(*t)[1], x0)
    return all_unit_distances(f42a_vertices(*t0)[0])

# F48A = genpetersen("f48a") but the resulting embedding is vertex-edge-degenerate, so...
def f48a():
    """Return a non-degenerate unit-distance embedding of the F48A graph."""
    R = (2 + 3*sqrt(2) + sqrt(12*sqrt(6)-26)) / 4
    r = (2 + 3*sqrt(2) - sqrt(12*sqrt(6)-26)) / 4
    L = R-1
    l = r-1
    u24 = unitroots(24)
    ring_R = [u*R for u in u24[::2]]
    ring_r = [u*r for u in u24[1::2]]
    ring_L = [u*L for u in u24[::2]]
    ring_l = [u*l for u in u24[1::2]]
    vertices = ring_R + ring_r + ring_L + ring_l
    edges = ring_edges(12, ((0, 1, 0), (0, 1, -1), (0, 2, 0), (1, 3, 0), (2, 3, 2), (2, 3, -3)))
    return (vertices, edges)

def f50a_vertices(t):
    u = root(1, 5, 1)
    table = {(): (0, 3, 4, 5, 6),
             (1,): (2, 35, 36, 37, 38),
             (1, 1, 1, 1): (48, 21, 22, 23, 24),
             (2,): (1, 28, 27, 26, 25),
             (1, 1): (34, 17, 18, 39, 40),
             (1, 2): (49, 46, 45, 44, 43),
             (1, 1, 1): (16, 19, 20, 41, 42),
             (2, 1): (33, 30, 29, 8, 7),
             (1, 1, 2): (47, 14, 13, 12, 11),
             (2, 1, 1): (15, 32, 31, 10, 9)}
    p0 = rect(star_radius(10), 0.9*pi)
    p3 = rect(star_radius(10, 3), -0.7*pi)
    p4 = cu(p3, -conj(p3))
    p5 = p4 + expj(t)
    p6 = cu(p5, u*p5)
    seeds = [p0, p3, p4, p5, p6]
    vertices = [None] * 50
    ops = {1: lambda z: u*z, 2: conj}
    for (aut, coset) in table.items():
        for (ring, i) in enumerate(coset):
            vertices[i] = -1j * reduce(lambda z, k: ops[k](z), aut, seeds[ring])
    return (vertices, re(vertices[40]) + 0.5)

def f50a():
    """Return a unit-distance embedding of the F50A graph, an embedding
    found by the computer (specifically the embedding_run() function in embeddingsearch)."""
    t0 = findroot(lambda t: f50a_vertices(t)[1], 2)
    return (f50a_vertices(t0)[0], lcf_edges(50, [21, -21, -19, 19, -19, 19, -19, 19, 21, -21]))

def f54a_vertices(t):
    u18 = unitroots(18)
    r0 = [u/2 for u in u18]
    z1 = cu(r0[1], r0[-1])
    r1 = [z1*u for u in u18]
    z2a = r1[0] + expj(t)
    z2b = circumcentre(z2a, u18[2]*z2a, r1[1])
    r2 = [u*z for u in unitroots(9) for z in (z2a, z2b)]
    vertices = r0 + r1 + r2
    return (vertices, abs(z2b - r1[1]) - 1)

def f54a():
    """Return a unit-distance embedding of the F54A graph."""
    t0 = findroot(lambda t: f54a_vertices(t)[1], 1.755)
    edges = ring_edges(18, ((0, 0, 9), (1, 0, 1), (1, 0, -1), (1, 2, 0), (2, 2, 1)))
    return (f54a_vertices(t0)[0], edges)

def f56a():
    """Return a unit-distance embedding of the F56A graph.
    Note that MathWorld's LCF notation for this is incorrect;
    it should be [11, 13, -13, -11]^14."""
    t = tan(pi/14)
    u = sqrt(polyval([-21, 98, 71], t*t))
    z1 = 2*sqrt(14*polyval([31*u, -20, -154*u, 104, 87*u, -68], t))
    z2 = 7*t*(t*t-3)**2 - 4*u
    a = (z1 + z2) / 64
    b = (z1 - z2) / 64

    u14 = unitroots(14)
    pa = mpc(a, 0.5)
    pb = mpc(b, 0.5)
    pac, pbc = conj(pa), conj(pb)
    d1 = abs(pa - u14[-1]*pb)**2 - 1
    d2 = abs(pb - u14[-2]*pa)**2 - 1
    vertices = [u*p for u in u14 for p in (pa, pb, pac, pbc)]
    return all_unit_distances(vertices)

def klein(a1=4.47, a2=2.42, a3=0.7, s1=1, s2=-1):
    """Return a unit-distance embedding of the cubic Klein graph (F56B)."""
    u7 = unitroots(7)
    z0 = star_radius(7)
    r0 = [z0*u for u in u7]
    z1 = z0 + expj(a1)
    z2 = z1 + expj(a2)
    z3 = z1 + expj(a3)
    r1 = [z1*u for u in u7]
    r2 = [z2*u for u in u7]
    r3 = [z3*u for u in u7]
    z4 = cu(*(r2[2], r3[0])[::s1])
    z5 = cu(*(r2[0], r3[1])[::s2])
    r4 = [z4*u for u in u7]
    r5 = [z5*u for u in u7]
    z6 = cu(0, r4[0], star_radius(7, 2), 1)
    z7 = cu(0, r5[0], star_radius(7, 3), 1)
    r6 = [z6*u for u in u7]
    r7 = [z7*u for u in u7]
    vertices = r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7
    edges = ring_edges(7, ((0, 0, 1), (0, 1, 0), (1, 2, 0), (1, 3, 0),
                           (2, 4, -2), (3, 4, 0), (2, 5, 0), (3, 5, -1),
                           (4, 6, 0), (5, 7, 0), (6, 6, 2), (7, 7, 3)))
    return (vertices, edges)

def f56c():
    """Return a unit-distance embedding of the F56C graph,
    the bipartite double cover of the Coxeter graph."""
    u14 = unitroots(14)
    z0 = star_radius(14, 5)
    r0 = [z0*u for u in u14]
    z1 = star_radius(14, 3)
    r1 = [z1*u for u in u14]
    z2 = cu(r1[4], r0[0])
    r2 = [z2*u for u in u14]
    z3 = cu(0, z2, star_radius(14), 1)
    r3 = [z3*u for u in u14]
    vertices = r0 + r1 + r2 + r3
    edges = ring_edges(14, ((0, 0, 5), (1, 1, 3), (2, 1, 4), (2, 0, 0), (2, 3, 0), (3, 3, 1)))
    return (vertices, edges)

def f60a(t=-0.35):
    """Return a unit-distance embedding of the F60A graph."""
    u15 = unitroots(15)
    z0 = star_radius(15, 7)
    r0 = [z0*u for u in u15]
    z1 = z0 + expj(t)
    r1 = [z1*u for u in u15]
    z2 = cu(r1[3], r1[0])
    r2 = [z2*u for u in u15]
    z3 = cu(0, z2, star_radius(15, 2), 1)
    r3 = [z3*u for u in u15]
    vertices = r0 + r1 + r2 + r3
    edges = ring_edges(15, ((0, 0, 7), (0, 1, 0), (2, 1, 0), (2, 1, 3), (2, 3, 0), (3, 3, 2)))
    return (vertices, edges)

def f62a_vertices(*params):
    u6 = unitroots(6)
    tree = [1j, 2j]
    tree.append(tree[-1] + expj(2.939))
    tree.append(tree[-1] + expj(-1.025))
    for (i, v) in enumerate((3, 1, 5, 6, 6)):
        tree.append(tree[v] + expj(params[i]))
    star = mpc(params[-2], params[-1])
    cc1 = circumcentre(tree[8],       star, u6[4]*tree[7])
    cc2 = circumcentre(u6[2]*tree[8], star,       tree[7])
    cc3 = circumcentre(u6[4]*tree[8], star, u6[2]*tree[7])
    cons = (abs(tree[2] - u6[1]*tree[5])**2 - 1,
            abs(tree[4] - u6[1]*tree[7])**2 - 1,
            abs(tree[3] - tree[8])**2 - 1,
            4*abs(tree[4])**2 - 1,
            abs(star - cc1)**2 - 1,
            abs(star - cc2)**2 - 1,
            abs(star - cc3)**2 - 1)
    vertices = [u*t for u in u6 for t in tree]
    vertices.extend(s*v for s in (1, -1) for v in (cc1, cc2, cc3, star))
    return (vertices, cons)

def f62a():
    """Return a unit-distance embedding of the F62A graph."""
    t0 = [-1.017, -0.819, 2.96, -0.282, -1.091, -0.624, 0.354]
    t0 = findroot(lambda *t: f62a_vertices(*t)[1], t0)
    return all_unit_distances(f62a_vertices(*t0)[0])

def f64a_vertices(a, b):
    u8 = unitroots(8)
    p1 = mpc(a, 0.5)
    p2 = mpc(b, 0.5)
    p3 = cu(u8[3]*p1, conj(p2), 2, 1)
    p4 = (u8[3]*p1 + p3) / 2
    d1 = abs(u8[1]*p3 - p4)**2 - 1
    d2 = abs(p1 - u8[1]*conj(p2))**2 - 1
    vertices = [u*p for u in u8 for p in (p1, p2, p3, p4)]
    vertices += [conj(p) for p in vertices]
    return vertices, (d1, d2)

def f64a():
    """Return a unit-distance embedding of the F64A graph."""
    t0 = findroot(lambda *t: f64a_vertices(*t)[1], (0.53, 1.6))
    return all_unit_distances(f64a_vertices(*t0)[0])

def f72a_vertices(t):
    u24 = unitroots(24)
    u12 = unitroots(12)
    z0 = star_radius(24, 11)
    r0 = [u*z0 for u in u24]
    z1 = star_radius(24, 7) * expj(t)
    r1 = [u*z1 for u in u24]
    z2 = cu(r0[0], r1[9])
    r2 = [u*z2 for u in u12]
    z3 = cu(r0[15], r1[6])
    r3 = [u*z3 for u in u12]
    vertices = r0 + r1 + r2 + r3
    return (vertices, abs(z2 - z3) - 1)

def f72a():
    """Return a unit-distance embedding of the F72A graph."""
    t0 = findroot(lambda t: f72a_vertices(t)[1], 2.2)
    return all_unit_distances(f72a_vertices(t0)[0])

def f74a_vertices(*params):
    u6 = unitroots(6)
    tree = [1j, 2j, 2j-expj(-pi/6), 2j+expj(pi/6)]
    params = [-1.04, 3.92] + list(params)
    for (i, v) in enumerate((2, 3, 4, 5, 5, 6, 7)):
        tree.append(tree[v] + expj(params[i]))
    star = mpc(params[-2], params[-1])
    cc1 = circumcentre(tree[8],       star,       -tree[9])
    cc2 = circumcentre(u6[2]*tree[8], star, -u6[2]*tree[9])
    cc3 = circumcentre(u6[4]*tree[8], star, -u6[4]*tree[9])
    cons = (abs(tree[6] - u6[1]*tree[8])**2 - 1,
            abs(tree[4] - tree[7])**2 - 1,
            abs(tree[9] + tree[10])**2 - 1,
            4*abs(tree[10])**2 - 1,
            abs(star - cc1)**2 - 1,
            abs(star - cc2)**2 - 1,
            abs(star - cc3)**2 - 1)
    vertices = [u*t for u in u6 for t in tree]
    vertices.extend(s*v for s in (1, -1) for v in (star, cc1, cc2, cc3))
    return (vertices, cons)

def f74a():
    """Return a unit-distance embedding of the F74A graph."""
    t0 = [2.91, 4.74, 5.5, 4.88, 5, -0.05, 0.07]
    t0 = findroot(lambda *t: f74a_vertices(*t)[1], t0)
    return all_unit_distances(f74a_vertices(*t0)[0])

def f78a_vertices(a, b, c):
    u13 = unitroots(13)
    pa = mpc(a, 0.5)
    pb = mpc(b, 0.5)
    pc = mpc(c, 0.5)
    pac, pbc, pcc = (conj(p) for p in (pa, pb, pc))
    d1 = abs(pa - u13[5]*pbc)**2 - 1
    d2 = abs(pb - u13[2]*pcc)**2 - 1
    d3 = abs(pc - u13[6]*pac)**2 - 1
    vertices = [u*p for u in u13 for p in (pa, pb, pc, pac, pbc, pcc)]
    return (vertices, (d1, d2, d3))

def f78a():
    """Return a unit-distance embedding of the F78A graph."""
    t0 = findroot(lambda *t: f78a_vertices(*t)[1], (-1.1, 1.6, 2.2))
    return all_unit_distances(f78a_vertices(*t0)[0])

def f80a(t=1.39):
    """Return a unit-distance embedding of the F80A graph."""
    u20 = unitroots(20)
    z0 = star_radius(20, 7)
    r0 = [z0*u for u in u20]
    z1 = z0 + expj(t)
    r1 = [z1*u for u in u20]
    z2 = cu(r1[2], r1[0])
    r2 = [z2*u for u in u20]
    z3 = cu(0, z2, star_radius(20, 3), 1)
    r3 = [z3*u for u in u20]
    vertices = r0 + r1 + r2 + r3
    edges = ring_edges(20, ((0, 0, 7), (0, 1, 0), (2, 1, 0), (2, 1, 2), (2, 3, 0), (3, 3, 3)))
    return (vertices, edges)

def f84a_vertices(p2, a, b, c):
    u7 = unitroots(7)
    p0 = star_radius(7)
    p1 = p0 + 1
    p3 = p2 + 1 # has a sign variation
    p4 = cu(p2, u7[3]*p1)
    p5 = mpc(a, 0.5)
    p6 = mpc(b, 0.5)
    p7 = mpc(c, 0.5)
    d1 = abs(p3 - u7[4]*p5)**2 - 1
    d2 = abs(p4 - u7[2]*p7)**2 - 1
    d3 = abs(p5 - u7[4]*conj(p6))**2 - 1
    d4 = abs(p6 - u7[-1]*p7)**2 - 1
    vertices = [u*p for u in u7 for p in (p0, p1, p2, p3, p4, p5, p6, p7)]
    vertices.extend([u*conj(p) for u in u7 for p in (p4, p5, p6, p7)])
    vertices = list(map(lambda z: z*1j, vertices))
    return (vertices, (d1, d2, d3, d4))

def f84a():
    """Return a unit-distance embedding of the F84A graph - not degenerate
    despite its looks. The graph is notable in having the simple PSL(2,8)
    as its automorphism group."""
    t0 = findroot(lambda *t: f84a_vertices(*t)[1], (-0.46, -1.44, 0.25, 0.75))
    return all_unit_distances(f84a_vertices(*t0)[0])

def f86a_vertices(*params):
    u6 = unitroots(6)
    tree = [1j, 2j, 2j-expj(-pi/6), 2j+expj(pi/6)]
    params = [5.24451, 5.34434, 5.00597] + list(params)
    for (i, v) in enumerate((2, 3, 4, 5, 5, 6, 7, 8, 9)):
        tree.append(tree[v] + expj(params[i]))
    star = mpc(params[-2], params[-1])
    cc1 = circumcentre(tree[10],       star,       tree[11])
    cc2 = circumcentre(u6[2]*tree[10], star, u6[2]*tree[11])
    cc3 = circumcentre(u6[4]*tree[10], star, u6[4]*tree[11])
    cons = (abs(tree[6] - u6[1]*tree[8])**2 - 1,
            abs(tree[4] - tree[7])**2 - 1,
            abs(tree[12] - tree[10])**2 - 1,
            4*abs(tree[9])**2 - 1,
            abs(tree[12] - u6[4]*tree[11])**2 - 1,
            abs(star - cc1)**2 - 1,
            abs(star - cc2)**2 - 1,
            abs(star - cc3)**2 - 1)
    vertices = [u*t for u in u6 for t in tree]
    vertices.extend(s*v for s in (1, -1) for v in (star, cc1, cc2, cc3))
    return (vertices, cons)

def f86a():
    """Return a unit-distance embedding of the F86A graph."""
    t0 = [3.60383, 3.44007, 4.34048, 5.63174, 3.26345, 0.488743, 0.113378, 0.236693]
    t0 = findroot(lambda *t: f86a_vertices(*t)[1], t0)
    return all_unit_distances(f86a_vertices(*t0)[0])

def foster_vertices(n, t):
    s2, s3 = (n&2)-1, ((n&1)<<1)-1
    c = star_radius(10)
    cp = c*root(1, 30, 1)
    a = rect(sec(pi/10), t)
    ap = rect(tan(pi/10), t+2*pi/5)
    b = cu(*(a, c)[::s2])
    bp = cu(*(ap, cp)[::s3])
    arc = (cp, bp, ap, a, b, c)
    vertices = [u*p for u in unitroots(15) for p in arc]
    return (vertices, abs(vertices[1] - vertices[82]) - 1)

def foster(i=5):
    """Return any one of six (depending on 0 <= i <= 5) unit-distance
    embeddings of the Foster graph (F90A)."""
    n, t0 = [(0, 0.38), (1, 1.35), (2, 0.15), (2, 1.18), (2, 4.68), (3, [1.5, 1.6])][i]
    tstar = findroot(lambda t: foster_vertices(n, t)[1], t0)
    return (foster_vertices(n, tstar)[0], lcf_edges(90, (17, -9, 37, -37, 9, -17)))

def foster_old_vertices(r):
    v3a = 0.265
    v3 = v3a * root(1, 5, 2)
    v2 = cu(v3, v3a)
    v5r = root(1, 20, 7) * r
    v5r2 = -v5r.conjugate()
    v5 = v5r * root(1, 15, 14)
    v0 = cu(v5r, v5r2)
    v1 = cu(v2, v0)
    v4 = cu(v3, v5)
    vgens = (v0, v1, v2, v3, v4, v5)
    vertices = [v*u for v in vgens for u in unitroots(15)]
    return (vertices, abs(v1 - v4*root(1, 15, 2)) - 1)

def foster_old():
    """Return the unit-distance embedding of the Foster graph (F90A)
    originally in Dounreay."""
    r0 = findroot(lambda r: foster_old_vertices(r)[1], 0.35)
    vertices = foster_old_vertices(r0)[0]
    edges = ring_edges(15, ((0, 1, 0), (1, 2, 0), (2, 3, 0), (3, 4, 0), (4, 5, 0), (5, 0, -1),
                            (0, 5, -2), (2, 3, -6), (4, 1, -2)))
    return (vertices, edges)

def f96a_vertices(a):
    u12 = unitroots(12)
    p1 = a+0.5j
    p2 = 2.3+0.5j
    p3 = cu(u12[-1]*p1, p2) # ring 3
    p4 = cu(p2, u12[1]*p1) # ring 2
    vertices = [u*p for u in u12 for p in (p1, p2, p3, p4)]
    vertices.extend([conj(p) for p in vertices])
    return (vertices, abs(p3 - u12[3]*conj(p4)) - 1)

def f96a():
    """Return a unit-distance embedding of the F96A graph."""
    a0 = findroot(lambda a: f96a_vertices(a)[1], 0.7)
    return all_unit_distances(f96a_vertices(a0)[0])

def f96b(a=2.32, b=1.92, c=-0.26, s1=-1, s2=1, s3=1, s4=-1):
    """Return a unit-distance embedding of the F96B graph."""
    u12 = unitroots(12)
    z2 = star_radius(12, 5)
    r2 = [u*z2 for u in u12]
    z3 = z2 + expj(a)
    r3 = [u*z3 for u in u12]
    z1 = z3 + expj(b)
    r1 = [u*z1 for u in u12]
    z4 = cu(*(u12[3]*z1, z3)[::s1])
    r4 = [u*z4 for u in u12]
    z5 = z1 + expj(c)
    r5 = [u*z5 for u in u12]
    z6 = cu(*(u12[-5]*z5, z4)[::s2])
    r6 = [u*z6 for u in u12]
    z8 = cu(*(u12[-4]*z6, z5)[::s3])
    r8 = [u*z8 for u in u12]
    z7 = cu(z8, 0, 1, star_radius(12)) if s4 == 1 else cu(0, z8, star_radius(12), 1)
    r7 = [u*z7 for u in u12]
    vertices = r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8
    edges = ring_edges(12, ((1, 1, 5), (1, 2, 0), (2, 0, 0), (3, 0, 3),
                            (3, 2, 0), (4, 0, 0), (5, 4, -5), (5, 3, 0),
                            (7, 5, -4), (7, 4, 0), (6, 7, 0), (6, 6, 1)))
    return (vertices, edges)

def f98a_vertices(a, t2, t3, t4, t5):
    u7 = unitroots(7)
    p1 = mpc(a, 0.5)
    p6 = 1.81+0.5j
    p7 = 2.18+0.5j
    p2 = p1 + expj(t2)
    p3 = p6 + expj(t3)
    p4 = p7 + expj(t4)
    p5 = p7 + expj(t5)
    cons = (abs(p1 - u7[3]*conj(p4))**2 - 1,
            abs(p2 - u7[-2]*conj(p3))**2 - 1,
            abs(p2 - u7[3]*p5)**2 - 1,
            abs(p3 - u7[-1]*p5)**2 - 1,
            abs(p4 - u7[1]*conj(p6))**2 - 1)
    vertices = [u*p for u in u7 for p in (p1, p2, p3, p4, p5, p6, p7)]
    vertices.extend([conj(p) for p in vertices])
    return (vertices, cons)

def f98a():
    """Return a unit-distance embedding of the F98A graph."""
    t0 = findroot(lambda *t: f98a_vertices(*t)[1], (-0.075, -2.3, -2.5, -2.8, -2.5))
    return all_unit_distances(f98a_vertices(*t0)[0])

@fixparams_unitdist(3.2, -2.5, 2.7, -2.7, -3.2)
def f98b(a, b, c, d, e):
    u = root(1, 7, 1)
    z1 = 3.175+0.5j
    z2 = z1+expj(a)
    z3 = z1+expj(b)
    z4 = z3+expj(c)
    z5 = z3+expj(d)
    z6 = z4+expj(e)
    z7 = -0.0375+0.5j
    cons = (abs(z2 - conj(z2)*u) - 1,
            abs(z2 - z5*u) - 1,
            abs(z4 - z7/u**2) - 1,
            abs(z5 - z6/u) - 1,
            abs(z6 - conj(z7)/u) - 1)
    vertices = symmetrise((z1, z2, z3, z4, z5, z6, z7), "D7")
    return (vertices, cons)

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
