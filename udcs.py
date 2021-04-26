#!/usr/bin/env python3
from shibuya import draw
from shibuya.cubesym import *
from shibuya.cubesym2 import *
from shibuya.igraph import genpetersen

for func in (heawood, pappus, f26a, coxeter, tutte8, dyck, f38a, f40a, f42a, f48a, f50a,
        f54a, f56a, klein, f56c, f60a, f62a, f64a, f72a, f74a, f78a, f80a, f84a, f86a,
        foster, f96a, f96b, f98a, f98b, biggssmith, f104a, f108a, f110a, f112a, f112b,
        f112c, f114a, f120a, f120b):
    draw(func(), func.__name__)

for s in ("cube", "petersen", "mobiuskantor", "dodecahedron", "desargues", "nauru"):
    draw(genpetersen(s), s)
