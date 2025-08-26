"""
Other unit-distance graphs...
"""
from mpmath import *
from shibuya.generators import (cu, star_radius, ring_edges,
        all_unit_distances, fixparams_unitdist, symmetrise)

def franklin():
    """Return a unit-distance embedding of the Franklin graph."""
    s2 = polyroots([9, -15, 16, -15, 9])[1]
    s3 = sqrt(3) * 1j
    s1 = cu(s2, s3)
    third = [v - s3/3 for v in (s3, s1, s2, 0)]
    vertices = [v*u for v in third for u in unitroots(3)]
    edges = ring_edges(3, ((0, 1, 0), (0, 3, 1), (0, 3, -1), (1, 2, 0), (1, 2, -1), (2, 3, 0)))
    return (vertices, edges)

@fixparams_unitdist()
def smallest_zerosym():
    """Return a unit-distance embedding of the smallest zero-symmetric graph
    (vertex-transitive, edges are partitioned into three orbits), which has LCF
    notation [5, -5]^9."""
    a = polyroots([16, -16*sqrt(3), 32, 12*sqrt(3), -1])[0]
    z1 = root(1, 12, 5)
    z2 = mpc(a, 0.5)
    z3 = cu(1j, z2)
    return symmetrise((z1, z2, z3), "D3")

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

@fixparams_unitdist()
def holt():
    """Return the unique non-degenerate unit-distance embedding of the Holt graph,
    the smallest half-transitive graph, with D9 symmetry."""
    def r(i):
        t = 4*cos(2**i*pi/9)**2
        return sqrt(polyroots([1, 1-t, t, -t/3])[0])
    return symmetrise((r(0), -r(1), -r(2)), "C9")

@fixparams_unitdist()
def gray():
    """Return a unit-distance embedding of the Gray graph, the smallest cubic
    edge- but not vertex-transitive graph with LCF [-25, 7, -7, 13, -13, 25]^9."""
    line = [p+d for p in polyroots([1, 0, -3, 1]) for d in (0, 1)]
    return symmetrise(line, "C9")

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

@fixparams_unitdist()
def balaban10():
    """Return a unit-distance embedding of the Balaban 10-cage."""
    p1 = 0.5j
    p2 = 0.8
    p3 = -abs(cu(root(1,5,2)*p2, p2))
    p4 = cu(p1, p2)
    p5 = cu(p3, p1)
    p6 = cu(p4, 0, 1, star_radius(10))
    p7 = cu(0, p5, star_radius(10), 1)
    return symmetrise((p1, p2, p3, p4, p5, p6, p7), "C10")
