# Draws a unit-distance embedding of the bidiakis cube.

house = (0, 1j, root(1, 6, 1), root(1, 6, 1) + 1j, 1, 1+1j)
house_edges = [(0, 1), (2, 3), (4, 5), (0, 2), (2, 4), (1, 3), (3, 5)]
left_house = [v * sqrt(1j) for v in house]
right_house = [(v - 1) * sqrt(-1j) + 1 for v in house]
right_house_edges = [(a + 6, b + 6) for (a, b) in house_edges]
vertices = left_house + right_house
edges = house_edges + right_house_edges + [(1, 6), (4, 11), (5, 7), (0, 10)]
