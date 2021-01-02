"""
Other unit-distance graphs...
"""
from mpmath import *
from shibuya.generators import cu, star_radius, ring_edges, all_unit_distances, delete_vertices

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

def blanusa2():
    """Draws a unit-distance embedding of the second Blanuša snark."""
    u4 = unitroots(4)
    v1 = [u/2 for u in u4]
    v2_s = cu(v1[1], v1[0])
    v2 = [v2_s*u for u in u4]
    vertices1 = v1 + v2
    edges1 = [(0, 2), (1, 3), (0, 4), (4, 1), (1, 5), (5, 2), (2, 6), (6, 3), (3, 7), (7, 0)]
    tvec = mpc(sqrt(111), 15) / sqrt(336)
    vertices2 = [v + tvec for v in vertices1]
    edges2 = [(a + 8, b + 8) for (a, b) in edges1]
    vertices = vertices1 + vertices2
    va = cu(vertices[4], vertices[12])
    vb = cu(vertices[13], vertices[5])
    vertices.extend([va, vb])
    edges = edges1 + edges2 + [(6, 14), (7, 15), (4, 16), (12, 16), (5, 17), (13, 17), (16, 17)]
    return (vertices, edges)

def franklin():
    """Return a unit-distance embedding of the Franklin graph."""
    s2 = polyroots([9, -15, 16, -15, 9])[1]
    s3 = sqrt(3) * 1j
    s1 = cu(s2, s3)
    third = [v - s3/3 for v in (s3, s1, s2, 0)]
    vertices = [v*u for v in third for u in unitroots(3)]
    edges = ring_edges(3, ((0, 1, 0), (0, 3, 1), (0, 3, -1), (1, 2, 0), (1, 2, -1), (2, 3, 0)))
    return (vertices, edges)

def harborth():
    """Return the canonical presentation of the Harborth graph. Based on Gerbracht (2007),
    Minimal Polynomials for the Coordinates of the Harborth Graph,
    https://arxiv.org/abs/math/0609360"""
    # Degree-22 polynomials for (what in Gerbracht's paper are) the y-coordinate of H and x-coordinate of C, but we turn the graph 90°
    Hsquare_pol = [427098112, -27098808320, 615643279360, -4779985142784, -32095868573376, 823044986987616, -3667116898760364, -27343071784237320, 273168911377174014, -441020584930952232, -123412000423046805, -12148787578527675]
    Csquare_pol = [427098112, -50083921920, 2572257472512, -77064294460416, 1563610131071808, -24051159678783648, 292733387369474292, -2623723693990622868, 15170804748275250138, -49933201015710366166, 83653148035178006805, -55268097000787592100]
    with extradps(100):
        H = sqrt(polyroots(Hsquare_pol)[2])
        C = 1j * sqrt(polyroots(Csquare_pol)[1])
    D = cu(C, H, 3, 2)
    B = C + (D-C)/3 * root(1, 6, 5)
    E = D + (C-D)/3 * root(1, 6, 1)
    G = cu(H, D, 2, 2)
    F = cu(G, E)
    J = chop(cu(G, F))
    A = chop(cu(F, B))
    aux_points = [(C+2*D)/3, (2*C+D)/3, (B+E)/2, (D+H)/2, (D+G)/2, (H+G)/2]
    quadrant = [A, B, C, D, E, F, G, H, J] + aux_points
    half = quadrant + [-p.conjugate() for p in quadrant if p.real > 0]
    vertices = half + [p.conjugate() for p in half if p.imag > 0]
    return all_unit_distances(vertices)

def rigid_heptagon(suppress=None):
    """Return a unit-distance graph containing a regular heptagon that was conjectured by
    Ed Pegg to be rigid in https://math.stackexchange.com/q/3954719/357390.
    I proved it to be rigid, even when two adjacent vertices are deleted;
    use (2,3) or (2,2) or (3,3) as the parameter to see these reduced rigid graphs."""
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

def mcgee():
    """Return a unit-distance embedding of the McGee graph.
    From https://math.stackexchange.com/q/1484002/357390"""
    u8 = unitroots(8)
    r0 = [u/2 for u in u8[::2]]
    r1 = [u/2 for u in u8[1::2]]
    pol = [32, -192, 624, -1344, 2496, -4896, 9664, -15360, 17838, -14220, 7425, -2268, 324]
    z0 = polyroots(pol)[-2] / (1-1j)
    z1 = cu(z0, r0[0])
    r2 = [z0*u for u in u8[::2]]
    r3 = [z1*u for u in u8[::2]]
    r4 = [z0.conjugate()*u for u in u8[::2]]
    r5 = [z1.conjugate()*u for u in u8[::2]]
    vertices = r0 + r1 + r2 + r3 + r4 + r5
    edges = ring_edges(4, ((0, 0, 2), (1, 1, 2), (1, 2, 3), (1, 4, 2), (0, 3, 0), (0, 5, 0),
                           (2, 3, 0), (2, 3, 1), (4, 5, 0), (4, 5, 3)))
    return (vertices, edges)

def khodulyov_pentagon():
    """Return Andrei Khodulyov's regular pentagon made rigid with 31 edges."""
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
