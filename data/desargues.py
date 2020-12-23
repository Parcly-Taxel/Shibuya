# Draws a unit-distance embedding of the Desargues graph.

outer = [u * subtend_r(1, 10) for u in unitroots(10)]
inner = [u * (subtend_r(1, 10) - 1) for u in unitroots(10)]
vertices = outer + inner
edges = [(i, (i + 1) % 10) for i in range(10)] + [(i + 10, (i + 3) % 10 + 10) for i in range(10)] + [(i, i + 10) for i in range(10)]
