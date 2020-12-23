# Draws the (edge-degenerate) unit-distance embedding of the Hamming (3,3) graph.
# Note that it has two degrees of freedom despite all the edges.

s1 = subtend_r(1, 9)
#r1 = [s1 * u for u in unitroots(9)]
s2 = subtend_r(2, 9)
#r2 = [s2 * u for u in unitroots(9)]
s3 = -subtend_r(4, 9)
#r3 = [-s3 * u for u in unitroots(9)]

pairs = {"1,2": (s1, s2), "1,3": (s1, s3), "2,3": (s2, s3)}

for (k, v) in pairs.items():
    a, b = v
    for (n, u) in enumerate(unitroots(9)):
        z = abs(a - u*b)
        print(k, n, findpoly(z, 6, maxcoeff=10000, tol=1e-20))

#vertices = r1 + r2 + r3
#edges = ring_edges(9, ((0, 0, 1), (1, 1, 2), (0, 1, 1), (0, 1, -1), (2, 2, 4), (2, 1, 2), (2, 1, -2), (0, 2, 4), (0, 2, -4)))
