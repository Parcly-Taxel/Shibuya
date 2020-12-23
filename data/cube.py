# Draws a unit-distance embedding of the cube skeleton.
square = [0, 1, 1+1j, 1j]
tvec = mpc(sqrt(0.5), sqrt(0.5))
square2 = [v + tvec for v in square]
vertices = square + square2
edges = ((0, 1), (1, 2), (2, 3), (3, 0),
         (4, 5), (5, 6), (6, 7), (7, 4),
         (0, 4), (1, 5), (2, 6), (3, 7))
