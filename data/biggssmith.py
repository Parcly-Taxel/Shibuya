# Draws a unit-distance embedding of the Biggs-Smith graph.

# Arbitrarily place the 1 and 4 rings up and the 8 ring down;
# these determine all other vertices

u17 = unitroots(17)
r1 = [subtend_r(1, 17) * u * 1j for u in unitroots(17)]
r4 = [subtend_r(4, 17) * u * 1j for u in unitroots(17)]
r8 = [subtend_r(8, 17) * u * -1j for u in unitroots(17)]
sh1 = cu(r1[0], r4[0])
rh1 = [sh1 * u for u in unitroots(17)]
sh2 = cu(sh1, r8[7])
rh2 = [sh2 * u for u in unitroots(17)]
s2 = cc(sh2, 0, 1, subtend_r(2, 17))
r2 = [s2 * u for u in unitroots(17)]

vertices = r1 + r4 + rh1 + r8 + rh2 + r2

edges = [(i, i + 34) for i in range(17)] + [(i + 17, i + 34) for i in range(17)] + \
        [(i, (i + 1) % 17) for i in range(17)] + [(i + 17, (i + 4) % 17 + 17) for i in range(17)] + \
        [(i + 51, (i + 8) % 17 + 51) for i in range(17)] + [(i + 85, (i + 2) % 17 + 85) for i in range(17)] + \
        [(i + 34, i + 68) for i in range(17)] + [(i + 68, (i + 7) % 17 + 51) for i in range(17)] + [(i + 68, i + 85) for i in range(17)]
