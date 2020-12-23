# Draws a unit-distance embedding of the Dyck graph.

r1 = unitroots(8)
r2 = [sqrt(2) * r for r in r1]
s3 = cc(r1[1], 0, 1, subtend_r(1, 8))
r3 = [s3 * r for r in r1]
s4 = cc(0, r2[0], subtend_r(3, 8), 1)
r4 = [s4 * r for r in r1]

vertices = r1 + r2 + r3 + r4
edges = [(i, (i + 1) % 8 + 8) for i in range(8)] + [(i, (i - 1) % 8 + 8) for i in range(8)] + \
        [(i, (i - 1) % 8 + 16) for i in range(8)] + [(i + 16, (i + 1) % 8 + 16) for i in range(8)] + \
        [(i + 8, i + 24) for i in range(8)] + [(i + 24, (i + 3) % 8 + 24) for i in range(8)]
