# Draws the smallest-largest-edge integral embedding (7) of the octahedral graph,
# from Harborth and Moller, Minimum Integral Drawings of the Platonic Graphs,
# Mathematics Magazine vol. 67 no. 5 (December 1994), 355-358.

v1 = mpf(13) / 3 * unitroots(6)[1] - mpf(2) / 3
v2 = mpc(-v1.real, v1.imag)
omega = unitroots(3)[1]
v3 = v1 * omega
v4 = v2 * omega
v5 = v1 / omega
v6 = v2 / omega
vertices = (v1, v2, v3, v4, v5, v6)

edges = ((0, 1), (0, 2), (0, 4), (0, 5), (1, 2), (1, 3), (1, 5), (2, 3), (2, 4), (3, 4), (3, 5), (4, 5))
uscale = 4
