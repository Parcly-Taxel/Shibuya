# Draws a unit-distance embedding of the Nauru graph.
# Embedding from Žitnik, Horvat and Pisanski,
# All Generalized Petersen Graphs are Unit-Distance Graphs, http://preprinti.imfm.si/PDF/01109.pdf

# Four radii, outer to inner
R = (6 + 7 * sqrt(3) + sqrt(15)) / 12
r = (6 + 7 * sqrt(3) - sqrt(15)) / 12
L = R - 1
l = r - 1
u12 = unitroots(12)
ring_R = [u * R for u in u12[::2]]
ring_r = [u * r for u in u12[1::2]]
ring_L = [u * L for u in u12[::2]]
ring_l = [u * l for u in u12[1::2]]
vertices = ring_R + ring_r + ring_L + ring_l

edges = [(i, i + 6) for i in range(6)] + [(i, (i - 1) % 6 + 6) for i in range(6)] + [(i, i + 12) for i in range(6)] + \
        [(i + 12, (i + 2) % 6 + 18) for i in range(6)] + [(i + 12, (i - 3) % 6 + 18) for i in range(6)] + [(i + 6, i + 18) for i in range(6)]
