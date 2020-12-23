# Draws a unit-distance embedding of the flower snark J_5.

s1 = 2j * subtend_r(2, 5)
r1 = [s1 * u for u in unitroots(5)]
s2 = (r1[1] + r1[4]) / 2
r2 = [s2 * u for u in unitroots(5)]
s3 = cu(s2, s1)
r3 = [s3 * u for u in unitroots(5)]
s4 = cc(0, s3, subtend_r(1, 5), 1)
r4 = [s4 * u for u in unitroots(5)]
vertices = r1 + r2 + r3 + r4
edges = ring_edges(5, ((2, 0, 0), (2, 3, 0), (2, 1, 0), (3, 3, 1), (1, 0, 1), (1, 0, -1)))
