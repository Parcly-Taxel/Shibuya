# Draws a unit-distance embedding of the Tutte 8-cage.
from itertools import combinations

def fingap(r2, vertices=False):
    r1 = 0.25
    u1 = rect(r2, 7 * pi / 10)
    u2 = rect(r2, 3 * pi / 10)
    v1 = cu(u1, u2) # ring 6
    
    u1 = rect(r1, 9 * pi / 10)
    u2 = rect(r1, -3 * pi / 10)
    v2 = cu(u1, u2) # ring 5
    
    u1 = rect(r1, -7 * pi / 10) # ring 1
    u2 = -r2 * 1j # ring 2
    v3 = cu(u1, u2) # ring 4
    
    v4 = cu(v2, v1) # ring 3
    
    if vertices: return (v1, v2, v3, v4, u1, u2)
    return abs(v3 - v4) - 1

fifth = fingap(findroot(fingap, 0.3), True)
vertices = [v * ur for v in fifth for ur in unitroots(5)]
edges = [p for p in combinations(range(30), 2) if almosteq(abs(vertices[p[0]] - vertices[p[1]]), 1)]
uscale = 0.5
