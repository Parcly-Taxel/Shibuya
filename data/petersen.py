# Draws a unit-distance embedding of the Petersen graph.

outer = [subtend_r(1, 5) * expjpi(mpf(2 * i) / 5) * 1j for i in range(5)]
inner = [-subtend_r(2, 5) * expjpi(mpf(2 * i) / 5) for i in range(5)]
vertices = outer + inner
edges = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
         (5, 7), (7, 9), (9, 6), (6, 8), (8, 5),
         (0, 5), (1, 6), (2, 7), (3, 8), (4, 9))
