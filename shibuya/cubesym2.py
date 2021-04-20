"""
Cubic symmetric graphs with more than 102 vertices.
"""
from mpmath import *
from shibuya.generators import fixparams_unitdist, symmetrise

@fixparams_unitdist(0.37, 1.47)
def f104a(a, b):
    p1 = mpc(a, 0.5)
    p2 = mpc(b, 0.5)
    d1 = abs(p1 - p2*root(1, 26, 4)) - 1
    d2 = abs(p2 - p1/root(1, 26, 1)) - 1
    return (symmetrise((p1, p2), "D26"), (d1, d2))
