# Draws a unit-distance embedding of the truncated cubical graph.

s1 = root(1, 12, 1)
s2 = s1 + 1
s3 = s1 + root(1, 6, 1)
s4 = cc(s1, 0, 1, subtend_r(1, 3))
r1 = [s1 * u for u in unitroots(6)]
r2 = [s2 * u for u in unitroots(6)]
r3 = [s3 * u for u in unitroots(6)]
r4 = [s4 * u for u in unitroots(6)]
vertices = r1 + r2 + r3 + r4
edges = ring_edges(6, ((0, 1, 0), (0, 2, 0), (1, 2, 0), (2, 1, 1), (3, 3, 2), (0, 3, 0)))
