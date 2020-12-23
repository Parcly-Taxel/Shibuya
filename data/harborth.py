# Draws the canonical presentation of the Harborth graph.
# Based on Gerbracht (2007), Minimal Polynomials for the Coordinates of the Harborth Graph, https://arxiv.org/abs/math/0609360
from itertools import combinations

# Degree-22 polynomials for (what in Gerbracht's paper are) the y-coordinate of H and x-coordinate of C, but we turn the graph 90°
Hsquare_pol = [427098112, -27098808320, 615643279360, -4779985142784, -32095868573376, 823044986987616, -3667116898760364, -27343071784237320, 273168911377174014, -441020584930952232, -123412000423046805, -12148787578527675]
H = sqrt(polyroots(Hsquare_pol)[2]) # ~3.621
Csquare_pol = [427098112, -50083921920, 2572257472512, -77064294460416, 1563610131071808, -24051159678783648, 292733387369474292, -2623723693990622868, 15170804748275250138, -49933201015710366166, 83653148035178006805, -55268097000787592100]
C = 1j * sqrt(polyroots(Csquare_pol)[1]) # ~2.980
# Rest of the major points
D = cc(C, H, 3, 2)
B = C + (D - C) / 3 * root(1, 6, 5)
E = D + (C - D) / 3 * root(1, 6, 1)
G = cc(H, D, 2, 2)
F = cu(G, E)
J = chop(cu(G, F))
A = chop(cu(F, B))
# Auxiliary points that make the graph true matchstick
aux_points = [(C + 2 * D) / 3, (2 * C + D) / 3, (B + E) / 2, (D + H) / 2, (D + G) / 2, (H + G) / 2]
quadrant = [A, B, C, D, E, F, G, H, J] + aux_points
half = quadrant + [-p.conjugate() for p in quadrant if p.real > 0]
vertices = half + [p.conjugate() for p in half if p.imag > 0]
edges = [p for p in combinations(range(52), 2) if almosteq(abs(vertices[p[0]] - vertices[p[1]]), 1, abs_eps=1e-30)]
