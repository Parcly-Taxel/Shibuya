"""
Cubic symmetric graphs. Most of the embeddings realised here were taken from MathWorld.
"""
from mpmath import *
from functools import reduce
from shibuya.generators import cu, star_radius, ring_edges, lcf_edges
from shibuya.generators import all_unit_distances, circumcentre

# F4A = tetrahedron() or complete(4) (not unit-distance)
# F6A = circulant(6, (1, 3)) or mobiusladder(3) (not unit-distance)
# F8A = genpetersen("cube")
# F10A = genpetersen("petersen")

def heawood():
    """Return the symmetric unit-distance embedding of the Heawood graph (F14A)
    tucked away in Mathematica's GraphData."""
    P = [10485760, 78643200, 263192576, 543686656, 812777472, 942080000, 843317248, 552468480, 208879616, -31170560, -99213312, -76779520, -32795648, 7878144, 17269760, 16256512, 11392032, 4836080, 3014064, 361320, 69498, -165789]
    # v0 is the only real root of the above polynomial
    v0 = polyroots(P, maxsteps=1000)[0]
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
    p8 = cu(p4, p2) # 1
    p9 = cu(p5, p3) # 1
    p10 = cu(-p8, p7) # 1
    V = (A, B, C, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
    return [v*s for v in V for s in (1, -1)]

def f26a():
    """Return a unit-distance embedding of the F26A graph."""
    L = [225763882172416, 677291646517248, -168428484689920, 20022957145325568, 110425351460487168, -62271608227627008, -1121399379324829696, 109067800210833408, 7540122093258539008, -7773427560423096320, -68122575204409933824, 52458094831308111872, 538023325894893109248, 78711684450812166144, -2822235439482864205824, -2995029967245104644096, 8042320072275722354688, 17883673052584642527232, -4091384557677847359488, -44565059159935463940096, -34337101102504632257536, 36701690560501134167040, 70580403666521815358208, 14232955138192862263296, -47035732072645104214528, -35808690608913191615616, 5894307940195308883552, 16637521343183451078624, 4349659847836042273980, -2275416792867611843748, -1127656410109393840797, 66655446601201742988, 72730432421368293597]
    with extradps(10):
        t = findroot(lambda x: polyval(L, x), 0.4) / 2
    return all_unit_distances(f26a_vertices(t))

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
    """Return a unit-distance embedding of the Tutte 8-cage (F30A)."""
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

# TODO F62A

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

# TODO F74A

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

def heawood_gerbracht():
    """Return a unit-distance embedding of the Heawood graph.
    This is the first construction given in Gerbracht (2008),
    Eleven Unit Distance Embeddings of the Heawood Graph,
    https://arxiv.org/abs/0912.5395"""
    p5, l5, p7, l7, p2, l3 = 0, 1, 1+1j, 1+2j, 2j, 1j
    # Degree-79 polynomial for the x-coordinate of l4
    P = [82521703002365615643033600000, 152135800369825007098920960000, -2120259444356145889512456192000, -8175821639408563679884718899200, 11025799477301561380923592949760, 149189048927171391219263917572096, 341989984727973884867396338188288, -763800345871643605733535512788992, -6892489761595983453459595854256128, -14303114368662112977785429692643328, 31343179682405215504161837658819584, 254616663098419271111012531383618560, 477056183905245971917488031692938304, -898635822877066299154282762314323520, -6556557400356413683063078157405200320, -11463391775661584618715895715025904128, 14201705397119143149709337683063717104, 109385892925207478360122518287948266224, 200727376265061817580032667984094835280, -60617631026953339799378305296984824656, -1097279690260575519531876572540059803892, -2435243231716218466580115477132980137292, -1650827959998751884275421145646879272940, 4733784662326174469816987234959768253776, 17321709733106215547946139735151891780269, 28955348159492426037443729536713509636773, 21867253654523569285667250014704999794577, -30934416501269415569285918492882277029311, -152272756904971138353148344210050406803617, -325157218431048323421805399697113970403121, -451645674349824891937650532097542435080453, -343872926425618669220741202688368202345065, 286935408276107233753158822122577885606822, 1923952833473734147443634652898764867278198, 5180867575272248126071836905848828341927070, 9840137643451726574992603743314811193317886, 12720991312674958659494390034544200285598942, 6140751881298455069763046326781936407849238, -19570606574427556470966233073766236628787234, -67869160289243415287139367956058055810404822, -126213399593210124294769323126921742027688497, -164022007275895644197052849670790737036540873, -146045321267662575006252144965793560225509061, -57662664820854923809493690824194000968797973, 77130998985650655864689962089382720577858101, 213238754173051016042819729417269617854966165, 355269696471069385886754716351566237266830009, 664103288660372783854070699409333594554864741, 1533983381070251025995839971747580678500964852, 3550298683130683434462662943215234037891507412, 7274584518541872070335933586256322019748139172, 12902691890291653798206974719870993995735753540, 19954407479150801176386566760213350973570196080, 27201188778043412156622512508710379868716241416, 32963773017875955980864755706102737727974961688, 35716781564909427260214641236767872783162088204, 34722813139066795200081139797717329223025992699, 30346538554876120431728853077314517314609386819, 23863989284324858347511498529857889181323950967, 16888459659695355326863471817880692818622047623, 10751995223766688842173817330681019107518783545, 6152912915312070842952691100801887803370907305, 3160933625571584072448347845721693351127774301, 1455265549140319863369871581645012065857955441, 599083193386406195758633497190777431543886358, 219897164806211674807756610736580167553542758, 71715126275516155874072784490774971066237326, 20690863770430719393270631202371992513434414, 5253121604626527413065008160498494678879110, 1166012291532956694933924468283307736346382, 224472408717775611491021156521892619843158, 37109973679879574898320679050599920287450, 5203227805425306398124203036880713293545, 608930205226991194133708856923335926849, 58239553681851019741523172701651095197, 4422420653730204080254904433581059629, 255652807673380729611728470237761555, 10528063279784456967398200502468691, 273675328487397647237991825000783, 3348011046054687446588586894387]
    with extradps(100):
        xl4 = findroot(lambda x: polyval(P, x), [-0.74, -0.72], maxsteps=100)
    yl4 = sqrt(4 - (1-xl4)**2)
    l4 = mpc(xl4, yl4)
    p4 = (l4 + l5) / 2
    p6 = cu(l4, l7)
    p3 = cu(l3, l4)
    l2 = cu(p4, p2)
    l1 = p3 + 1
    l6 = cu(p5, p6)
    p1 = cu(l6, l1)
    vertices = (l1, l2, l3, l4, l5, l6, l7, p1, p2, p3, p4, p5, p6, p7)
    edges = ring_edges(7, ((0, 1, 0), (0, 1, -1), (0, 1, 2)))
    return (vertices, edges)
