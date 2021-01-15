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
    return (vertices, edges)

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

def mcgee(mode=0):
    """Return a unit-distance embedding of the McGee graph,
    from https://math.stackexchange.com/q/1484002/357390. mode (0 or 1)
    selects between two algebraically related forms."""
    a = polyroots([1, 0, -129, -218])[mode]
    t = polyroots([1, (a**2+3*a-52)/68, (a+3)/16])[1]
    core = [u/2 for u in unitroots(8)]
    z0 = mpc(0.5+t, -root(1-t**2, 2, mode))
    z1 = cu(z0, core[1])
    outer = [z*u for z in (z0, conj(z0), z1, conj(z1)) for u in unitroots(4)]
    vertices = core + outer
    return all_unit_distances(vertices)
