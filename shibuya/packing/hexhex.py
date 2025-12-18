"""
Packings of unit-side hexagons into the smallest possible hexagon.
Cf. https://erich-friedman.github.io/packing/hexinhex
and https://arxiv.org/abs/2511.02864 for the 12-hexagon packing.
"""
from mpmath import *
from shibuya.draw import drawing

def p12():
    w = unitroots(6) # rotations by multples of 60 degrees
    with extradps(100):
        t = sqrt(polyroots([27, -27, 63, -295, 276, -16])[0])
        u = (sqrt(4 - 3*t*t) - t)/2
        a = 2 / (t+u)
        s = sqrt(polyroots([27, -108, -8118, -211591, 6728352, -41165056])[0])
    h1 = (t*w[2], u)
    h2 = (t*w[2]+a, u+a)
    h3 = (t*w[2]+a*w[2], u+a*w[2])
    h4 = (t*w[2]+a*(1+w[2]), u+a*(1+w[2]))
    third = [(v1-s*w[1], v2-s*w[1]) for (v1, v2) in (h1, h2, h3, h4)]
    third1 = [(v1*w[2], v2*w[2]) for (v1, v2) in third]
    third2 = [(v1*w[4], v2*w[4]) for (v1, v2) in third]
    return (s, third + third1 + third2)

def hexpath(v1, v2):
    w = unitroots(6)
    v3 = v2 + w[1]*(v2-v1)
    v4 = v1 + 2*w[1]*(v2-v1)
    v5 = v2 + 2*w[2]*(v2-v1)
    v6 = v1 + w[2]*(v2-v1)
    l = [c for v in [v1, v2, v3, v4, v5, v6] for c in (v.real, v.imag)]
    return [["M"] + l, ["Z"]]

def draw_packing(data, outfn, scale=400):
    s, hexsides = data
    res = drawing(scale, (2*s, sqrt(3)*s), (-s, -sqrt(3)/2*s))
    for side in hexsides:
        res.add_path(hexpath(*side), {"fill": "#6dc6fb", "stroke": "#1c92cd", "stroke-width": 0.01})
    res.add_path(hexpath(s, s*root(1,6,1)), {"fill": "none", "stroke": "#000", "stroke-width": 0.01})
    res.write(outfn)
