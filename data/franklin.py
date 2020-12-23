# Draws a unit-distance embedding of the Franklin graph.

s2 = mpc((5 + sqrt(33)) / 12, (sqrt(75) - sqrt(11)) / 12)
s1 = cu(s2, sqrt(-3))
third = [v - sqrt(-3) / 3 for v in (sqrt(-3), s1, s2, 0)]
vertices = [u * v for u in unitroots(3) for v in third]
edges = [(i, (i + 1) % 12) for i in range(12)] + [(i, (i - 5) % 12) for i in range(0, 12, 4)] + [(i, (i - 3) % 12) for i in range(1, 12, 4)]
