# Draws the smallest-largest-edge integral embedding of the icosahedron skeleton,
# from Harborth and Moller, Minimum Integral Drawings of the Platonic Graphs,
# Mathematics Magazine vol. 67 no. 5 (December 1994), 355-358.

v0 = 0
v1 = 8
v2 = cc(v0, v1, 4, 8)
v3 = v2 + 8
v4 = cc(v0, v1, 7, 6)
v5 = v4 + 8
v7 = cc(v2, v3, 4, 6)
v6 = v7 - 8
v9 = cc(v4, v5, 4, 6)
v8 = v9 - 8
v10 = cc(v6, v7, 7, 6)
v11 = v10 + 8

vertices = (v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11)
edges = [(0, 1), (0, 2), (0, 4), (0, 6), (0, 8), (1, 2), (1, 3), (1, 4), (1, 5), (2, 3), (2, 6), (2, 7),
         (3, 5), (3, 7), (3, 11), (4, 5), (4, 8), (4, 9), (5, 9), (5, 11), (6, 7), (6, 8), (6, 10),
         (7, 10), (7, 11), (8, 9), (8, 10), (9, 10), (9, 11), (10, 11)]
uscale = 4
