# Draws a unit-distance embedding of the dodecahedron skeleton.

outer = [u * subtend_r(1, 10) for u in unitroots(10)]
inner_s = cc(0, outer[0], subtend_r(1, 5), 1)
inner = [inner_s * u for u in unitroots(10)]
vertices = outer + inner
edges = [(i, (i + 1) % 10) for i in range(10)] + [(i + 10, (i + 2) % 10 + 10) for i in range(10)] + [(i, i + 10) for i in range(10)]
