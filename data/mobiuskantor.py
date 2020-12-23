# Draws a unit-distance embedding of the Möbius–Kantor graph.

u16 = unitroots(16, primitive=True) # Conveniently gives all the angles, for both circles
vertices = [u * subtend_r(1, 8) for u in u16] + [u * subtend_r(3, 8) for u in u16]
edges = [(i, (i + 1) % 8) for i in range(8)] + [(i + 8, (i + 3) % 8 + 8) for i in range(8)] + [((i + 1) % 8, i + 8) for i in range(8)]
