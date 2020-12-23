# Draws a unit-distance embedding of the Coxeter graph.

u14 = unitroots(14)
ring1 = [u * subtend_r(3, 7) for u in u14[::2]]
ring2 = [u * subtend_r(2, 7) for u in u14[1::2]]
ring3_s = cu(ring2[0], ring1[0])
ring3 = [u * ring3_s for u in u14[::2]]
ring4_s = cc(0, ring3_s, subtend_r(1, 7), 1)
ring4 = [u * ring4_s for u in u14[::2]]
vertices = ring1 + ring2 + ring3 + ring4

edges = [(i, (i + 3) % 7) for i in range(7)] + [(i + 7, (i + 2) % 7 + 7) for i in range(7)] + [(i + 21, (i + 1) % 7 + 21) for i in range(7)] + \
        [(i, i + 14) for i in range(7)] + [(i + 7, i + 14) for i in range(7)] + [(i + 14, i + 21) for i in range(7)]
