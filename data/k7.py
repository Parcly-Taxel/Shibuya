# Draws the minimal integral embedding of K_7 (maximum edge length = 33, see https://oeis.org/A096872).
# K_8 and K_9 require a max-length of 56.

v0, v1 = -8, 8
v2 = cc(v0, v1, 32, 32)
v3 = cc(v0, v1, 16, 28)
v5 = cc(v0, v1, 24, 32)
v4 = -v3.conjugate()
v6 = -v5.conjugate()
vertices = (v0, v1, v2, v3, v4, v5, v6)
edges = [(a, b) for a in range(6) for b in range(a + 1, 7)]
uscale = 16
