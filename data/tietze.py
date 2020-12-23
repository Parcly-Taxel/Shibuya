# Draws a unit-distance embedding of Tietze's graph.

t = (13 * sqrt(3) + sqrt(13 * (40 * sqrt(3) - 9))) / 52 # see tietze-exact.py and tietze.singular for how this was derived
b = 1 / sqrt(3)
c = b + root(1, 3, 1)
c_ = c * root(1, 3, 1)
d = c - mpc(t, sqrt(1 - t * t))
e = cu(d, c_)
vertices = [v * u for u in unitroots(3) for v in (b, c, d, e)]
edges = [(0, 1), (1, 2), (2, 3), (0, 4), (3, 5), (2, 7),
         (4, 5), (5, 6), (6, 7), (4, 8), (7, 9), (6, 11),
         (8, 9), (9, 10), (10, 11), (8, 0), (11, 1), (10, 3)]
