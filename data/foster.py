# Draws a unit-distance embedding of the Foster graph.

def fingap(r, vertices=False):
    # Ring 3 radius = 0.265 looks aesthetically best
    v3_ = 0.265
    v3 = v3_ * root(1, 5, 2)
    v2 = cu(v3, v3_)
    
    v5r = root(1, 20, 7) * r
    v5r_ = -v5r.conjugate()
    v5 = v5r * root(1, 15, 14)
    v0 = cu(v5r, v5r_)
    
    v1 = cu(v2, v0)
    v4 = cu(v3, v5)
    
    if vertices: return (v0, v1, v2, v3, v4, v5)
    return abs(v1 - v4 * root(1, 15, 2)) - 1

fifteenth = fingap(findroot(fingap, 0.35), True)
vertices = [v / ur for ur in unitroots(15) for v in fifteenth]
edges = [(i, (i + 1) % 90) for i in range(90)] + \
        [(i, (i + 17) % 90) for i in range(0, 90, 6)] + \
        [(i, (i + 37) % 90) for i in range(2, 90, 6)] + \
        [(i, (i +  9) % 90) for i in range(4, 90, 6)]
uscale = 0.5
