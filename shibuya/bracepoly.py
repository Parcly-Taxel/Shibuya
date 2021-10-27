"""
Braced polygons, triangle-free and otherwise.
All graphs are Laman unless otherwise stated; if such a graph has v vertices
it has 2v-3 edges.
"""
from mpmath import *
from shibuya.generators import cu, star_radius, ring_edges, all_unit_distances, remove_edges, delete_vertices
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

def rigid_hexagon():
    """Return the trivial braced regular hexagon."""
    vertices = unitroots(6) + [0]
    return (vertices, ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
                       (6, 1), (6, 2), (6, 3), (6, 4), (6, 5)))

def rigid_heptagon(suppress=(2,3)):
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

def khodulyov_octagon():
    """Return Andrei Khodulyov's 17-vertex rigid regular octagon."""
    core = khodulyov_square()[0]
    z0 = core[0] + 1j
    core.extend([z0, cu(core[4], z0), cu(z0, core[4])])
    core.append(cu(core[1], z0))
    core.append(cu(core[4], core[-1]))
    core.append(cu(core[5], core[-1]))
    return all_unit_distances(core)

def pl_cell(a, b, c, d):
    """If in the chain abcd the angles abc and bcd are equal, return the extra
    vertices of a Peaucellier–Lipkin linkage enforcing that equality."""
    z0 = b+d-c
    z1 = a+d-c
    y0 = cu(z1, c, 1, sqrt(3))
    y1 = cu(c, z1, sqrt(3), 1)
    m = y0+y1-z1
    a0 = cu(y0, c)
    a1 = cu(c, y0)
    a2 = cu(y1, c)
    a3 = cu(c, y1)
    return [z0, z1, y0, y1, m, a0, a1, a2, a3]

def khodulyov_nonagon():
    """Return Andrei Khodulyov's 27-vertex rigid regular nonagon."""
    z0 = star_radius(9)
    outer = [z0*u for u in unitroots(9)]
    cell1 = pl_cell(outer[7], outer[8], outer[0], outer[1])
    cell2 = pl_cell(outer[2], outer[1], outer[0], outer[8])
    z1 = cu(outer[3], outer[5])
    vertices = outer + cell1 + cell2[1:] + [z1]
    return all_unit_distances(vertices)

@remove_edges(lambda e: e == (0, 18))
def khodulyov_decagon():
    """Return Andrei Khodulyov's 29-vertex rigid regular decagon."""
    core = khodulyov_pentagon()[0]
    origin = core[10]
    l0, h0 = core[-1], cu(core[7], core[1])
    l1, h1 = cu(origin, h0), cu(l0, h0)
    l2, h2 = cu(origin, h1), cu(l1, h1)
    l3, h3 = cu(origin, h2), cu(l2, h2)
    h4 = cu(l3, h3)
    l4 = cu(origin, h4)
    f0 = cu(core[9], l4)
    f1 = cu(f0, h4)
    f2 = cu(core[3], f0)
    vertices = core + (h0, l1, h1, l2, h2, l3, h3, l4, h4, f0, f1, f2)
    return all_unit_distances(vertices)

def rigid_hendecagon():
    """Return a unit-distance braced regular hendecagon (11-gon) with 41 vertices
    and 79 edges, a large improvement over Khodulyov's 155-edge bracing.
    This is based on the cyclotomic field decomposition
    (sqrt(-11)-1) / 2 = z + z^3 + z^4 + z^5 + z^9 with z the first primitive
    11th root of unity.
    The proof can be found in brace11gonproof.py and roughly follows my
    braced heptagon proof, except that the Jacobian is deficient by two ranks
    and the nullspace function's coordinates all have saddle points – the directions
    in which those saddles stay at zero are different, however, which establishes
    second-order rigidity."""
    z0 = star_radius(11)
    outer = [z0*u for u in unitroots(11)]
    z1i = outer[1] + outer[-1] - outer[0]
    z2i = outer[2] + outer[-1] - outer[0]
    z1 = z1i * root(1, 11, 7)
    z2 = z2i * root(1, 11, 7)
    z3 = outer[3] + outer[-1] - outer[0]
    z13 = cu(z1, z3)
    z31 = cu(z3, z1)
    z23 = cu(z2, z3)
    z32 = cu(z3, z2)
    spindle = [z1, z2, z3, z13, z31, z23, z32]
    spindles = [v*root(1, 11, -4*k) for k in range(4) for v in spindle]
    vertices = outer + [z1i, z2i] + spindles
    return all_unit_distances(vertices)

@remove_edges(lambda e: e == (14, 15) or e == (15, 25))
def khodulyov_dodecagon():
    """Return Andrei Khodulyov's 26-vertex rigid regular dodecagon."""
    core = khodulyov_square()[0]
    p0 = cu(core[1], core[0])
    p1 = cu(core[1], p0)
    p2 = cu(p0, core[0])
    p3 = cu(p1, p0)
    p4 = cu(p0, p2)
    p5 = cu(core[2], core[1])
    p6 = cu(core[0], core[3])
    p7 = cu(p5, p1)
    p8 = cu(p2, p6)
    p9 = cu(p7, p1)
    p10 = cu(p2, p8)
    p11 = cu(p9, p3)
    p12 = cu(p4, p10)
    p13 = cu(p11, p3)
    p14 = cu(p13, p12)
    vertices = core + [p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14]
    return all_unit_distances(vertices)

@remove_edges(lambda e: e == (20, 34))
def rigid_tridecagon():
    """Return a unit-distance braced regular tridecagon (13-gon) with 77 vertices
    and 151 edges. This can be proved rigid in an analogous manner to the 11-gon case –
    the nullspace here is 6-dimensional, but one coordinate of the nullspace function
    has a positive-definite Hessian matrix, from which the result follows."""
    p = [star_radius(13)*u for u in unitroots(13)]
    def bridge(i, d):
        q0 = p[i] + p[(i+2*d)%13] - p[(i+d)%13]
        q1 = p[i] + p[(i-2*d)%13] - p[(i-d)%13]
        q2 = q0 + q1 - p[i]
        q3 = p[(i-3*d)%13] + q1 - p[(i-2*d)%13]
        q4 = q2 + q3 - q1
        q5 = p[(i+3*d)%13] + p[(i+5*d)%13] - p[(i+4*d)%13]
        q6 = q5 + p[(i+6*d)%13] - p[(i+5*d)%13]
        rl0 = cu(*[q3, q4][::d])
        rr0 = cu(*[q4, q3][::d])
        rl1 = (q6 + rl0) / 2
        rr1 = (q6 + rr0) / 2
        rl2 = cu(*[q6, rl1][::d])
        rl3 = cu(*[rl1, rl0][::d])
        rr2 = cu(*[rr0, rr1][::d])
        rr3 = cu(*[rr1, q6][::d])
        return [[q0, q1, q2, q3, q4, q5, q6], [rl0, rl1, rl2, rl3, rr0, rr1, rr2, rr3]]
    bridges = [bridge(i, d) for (i, d) in ((5, 1), (5, -1), (0, 1), (0, -1), (-5, 1))]
    qs, rs = zip(*bridges)
    unique_qs = [qs[0][0]]
    for qbridge in qs:
        for q in qbridge:
            if min(abs(q-uq) for uq in unique_qs) > 1e-12:
                unique_qs.append(q)
    vertices = p + unique_qs + [r for rset in rs for r in rset]
    return all_unit_distances(vertices)

def rigid_tetradecagon():
    """Return Somsky's braced 14-gon (47 vertices, 91 edges), based on the rigid heptagon."""
    core = rigid_heptagon((2, 3))[0]
    core = [v - core[0] for v in core]
    z3, z2, z1 = core[6:3:-1]
    r3 = [z3*root(1, 14, n+1) for n in range(7)]
    r2 = [z2*root(1, 14, n+1) for n in range(9)]
    r1 = [z1*root(1, 14, n+1) for n in range(12)]
    vertices = core + r3 + r2 + r1
    return all_unit_distances(vertices)

def rigid_hexadecagon():
    """Return a braced regular 16-gon with 56 vertices and 109 edges."""
    rhombs = [0, 1, root(1,16,7), 1+root(1,16,1), root(1,16,7)+root(1,16,6),
              1+root(1,16,1)+root(1,16,2), root(1,16,7)+root(1,16,6)+root(1,16,5),
              1+root(1,16,1)+root(1,16,2)+root(1,16,3),
              root(1,16,7)+root(1,16,6)+root(1,16,5)+root(1,16,4)]
    for (i, j, k) in ((0,3,1), (0,4,2), (5,9,3), (6,10,4), (9,10,0), (12,13,10), (11,13,9),
                      (6,14,12), (5,15,11), (8,16,6), (7,17,5), (14,18,16), (15,19,17),
                      (13,20,14), (13,21,15), (22,23,13), (21,24,23), (20,24,22), (25,26,24),
                      (18,26,20), (19,25,21), (8,28,18), (7,29,19), (27,28,26), (27,29,25),
                      (30,32,28), (31,33,29), (32,33,27)):
        rhombs.append(rhombs[i] + rhombs[j] - rhombs[k])
    for (i, j) in ((3,15), (4,14), (9,17), (10,16), (18,24), (19,24)):
        rhombs.append(cu(rhombs[i], rhombs[j]))
        rhombs.append(cu(rhombs[j], rhombs[i]))
    rhombs.extend(pl_cell(*(rhombs[i] for i in (6, 4, 2, 0)))[2:])
    return all_unit_distances(rhombs)

@remove_edges(lambda e: 21 <= e[0] < e[1] <= 24 or 36 <= e[0] < e[1] <= 44)
def rigid_octadecagon():
    """Return a braced regular 18-gon with 60 vertices and 117 edges."""
    rhombs = [[0], [root(1, 18, 4), root(1, 18, -4)]]
    for k in range(3, -4, -1):
        middle = [rhombs[-1][i] + rhombs[-1][i+1] - rhombs[-2][i] for i in range(len(rhombs[-2]))]
        rhombs.append([rhombs[-1][0] + root(1, 18, k)] + middle + [rhombs[-1][-1] + root(1, 18, -k)])
    rhombs = [v for layer in rhombs+[[2*star_radius(18)]] for v in layer]
    rhombs.extend(pl_cell(rhombs[3], rhombs[1], rhombs[0], rhombs[2])[2:])
    rhombs.extend(pl_cell(rhombs[5], rhombs[2], rhombs[0], rhombs[1])[2:])
    return all_unit_distances(rhombs)

def khodulyov_polygon(n):
    """For n >= 7 construct the rigid regular n-gon through Khodulyov's
    equal-angle construction, which uses 19(n-3) + 4 - (n mod 2) edges."""
    z0 = star_radius(n)
    outer = [z0*u for u in unitroots(n)]
    vertices = outer[:]
    for k in range(2-(n//2), n//2, 2):
        cell = pl_cell(outer[k-2], outer[k-1], outer[k], outer[k+1])
        if 2*k+4 < n:
            cell += pl_cell(outer[k+2], outer[k+1], outer[k], outer[k-1])[1:]
        vertices.extend(cell)
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
