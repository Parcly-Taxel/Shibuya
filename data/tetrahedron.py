# Draws the smallest-largest-edge integral embedding of the tetrahedral graph,
# from Harborth and Moller, Minimum Integral Drawings of the Platonic Graphs,
# Mathematics Magazine vol. 67 no. 5 (December 1994), 355-358.

v1 = 0
v2 = 4
v3 = cc(v1, v2, 2, 4)
v4 = cc(v1, v2, 4, 2)
vertices = (v1, v2, v3, v4)
edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

uscale = 4
