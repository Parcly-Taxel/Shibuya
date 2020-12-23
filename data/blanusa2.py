# Draws a unit-distance embedding of the second Blanuša snark.

u4 = unitroots(4)
v1 = [0.5 * u for u in u4]
v2_s = cu(v1[1], v1[0])
v2 = [v2_s * u for u in u4]
vertices1 = v1 + v2
edges1 = [(0, 2), (1, 3), (0, 4), (4, 1), (1, 5), (5, 2), (2, 6), (6, 3), (3, 7), (7, 0)]
tvec = mpc(sqrt(mpf(111) / 336), sqrt(mpf(225) / 336))
vertices2 = [v + tvec for v in vertices1]
edges2 = [(a + 8, b + 8) for (a, b) in edges1]
vertices = vertices1 + vertices2

va = cu(vertices[4], vertices[12])
vb = cu(vertices[13], vertices[5])
vertices.extend([va, vb])

print(abs(vertices[10] - vertices[17]))

edges = edges1 + edges2 + [(6, 14), (7, 15), (4, 16), (12, 16), (5, 17), (13, 17), (16, 17)]
