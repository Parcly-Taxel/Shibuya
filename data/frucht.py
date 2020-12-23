# Draws a unit-distance embedding of the Frucht graph.
from itertools import combinations

frame = [0, 1, 1j, 1+1j, root(1, 12, 5), 1 + root(1, 12, 1), 1j + root(1, 6, 1)]
t = polyroots([53248, -39936, 19328, 44832, 9604, -11604, -13899, 738, -36423, -19440, 1296])[2]
D = frame[4] + mpc(sqrt(1 - t * t), t)
E = cu(D, frame[6])
F = cc(D, frame[5], sqrt(3), 1)
M = (E + F) / 2
N = M + D - E
vertices = frame + [D, E, F, M, N]
edges = [p for p in combinations(range(12), 2) if almosteq(abs(vertices[p[0]] - vertices[p[1]]), 1, abs_eps=1e-30) and p not in ((2, 3), (7, 10))]
