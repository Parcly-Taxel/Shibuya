# Draws a unit-distance embedding (in the plane!) of the truncated tetrahedral graph.

r1 = [0.5, 0.5j, -0.5, -0.5j]
cpsi = (-3 + 3 * sqrt(3) + sqrt(4 + 6 * sqrt(3))) / 8 # the cosine of an angle
v2 = 0.5 + mpc(cpsi, sqrt(1 - cpsi ** 2))
r2 = [v2 * 2 * r for r in r1]
v3 = 0.5 + mpc(cpsi, sqrt(1 - cpsi ** 2)) * root(1, 6, 1)
r3 = [v3 * 2 * r for r in r1]

vertices = [v * (r3[3].conjugate() / abs(r3[3])) for v in r1 + r2 + r3]
edges = ((0, 2), (1, 3),
         (0, 4), (1, 5), (2, 6), (3, 7),
         (0, 8), (1, 9), (2, 10), (3, 11),
         (4, 8), (5, 9), (6, 10), (7, 11),
         (8, 5), (9, 6), (10, 7), (11, 4))
