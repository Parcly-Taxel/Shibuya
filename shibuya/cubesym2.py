"""
Cubic symmetric graphs with more than 102 vertices. The lettering for graphs of the
same order follows Conder's list (https://www.math.auckland.ac.nz/~conder/symmcubic2048list.txt),
which sorts by arc-transitivity, not Foster's list.
"""
from mpmath import *
from shibuya.generators import cu, star_radius
from shibuya.generators import fixparams_unitdist, symmetrise, remove_edges

@fixparams_unitdist(0.37, 1.47)
def f104a(a, b):
    p1 = mpc(a, 0.5)
    p2 = mpc(b, 0.5)
    d1 = abs(p1 - p2*root(1, 26, 4)) - 1
    d2 = abs(p2 - p1/root(1, 26, 1)) - 1
    return (symmetrise((p1, p2), "D26"), (d1, d2))

@remove_edges(lambda e: {e[0]//9, e[1]//9} < {0, 2, 4} and e[0]//9 != e[1]//9)
@fixparams_unitdist(0.18, 3, 3)
def f108a(a, b, c):
    z1 = -star_radius(9)
    z2 = z1 - 1
    z3 = -star_radius(9, 2)
    z4 = z3 + 1
    z5 = star_radius(9, 4)
    z6 = z5 + 1
    z7 = z2 + expj(a)
    z8 = z4 + expj(b)
    z9 = z6 + expj(c)
    d1 = abs(z7 - z8*root(1, 9, 1)) - 1
    d2 = abs(z8 - conj(z9)*root(1, 9, -2)) - 1
    d3 = abs(z9 - conj(z7)*root(1, 9, 5)) - 1
    vertices = symmetrise((z1, z2, z3, z4, z5, z6), "C9") + symmetrise((z7, z8, z9), "D9")
    return (vertices, (d1, d2, d3))

@fixparams_unitdist(-1.5, -1.4, -0.55, -0.25, 0.8)
def f110a(a, b, c, d, e):
    u = unitroots(11)
    p1 = mpc(a, 0.5)
    p2 = mpc(b, 0.5)
    p3 = mpc(c, 0.5)
    p4 = mpc(d, 0.5)
    p5 = mpc(e, 0.5)
    d1 = abs(p1 - conj(p3)*u[-1]) - 1
    d2 = abs(p3 - conj(p2)*u[-3]) - 1
    d3 = abs(p2 - p4*u[2]) - 1
    d4 = abs(p4 - p5*u[5]) - 1
    d5 = abs(p5 - conj(p1)*u[-4]) - 1
    vertices = symmetrise((p1, p2, p3, p4, p5), "D11")
    return (vertices, (d1, d2, d3, d4, d5))

@fixparams_unitdist(1.1, 1.7, 0.7, 1.75)
def f112a(a, b, c, d):
    u = unitroots(14)
    v1 = mpc(a, 0.5)
    v2 = mpc(b, 0.5)
    h1 = mpc(0.5, c)
    h2 = mpc(0.5, d)
    d1 = abs(v1 - h1*u[1]) - 1
    d2 = abs(v1 - conj(h2)*u[5]) - 1
    d3 = abs(v2 - conj(h1)*u[2]) - 1
    d4 = abs(v2 - h2*u[-1]) - 1
    vertices = symmetrise((v1, v2, h1, h2), "D14")
    return (vertices, (d1, d2, d3, d4))

@fixparams_unitdist(1.5, -1.4, 0.3, 2.5, 2.7)
def f112b(a, b, c, d, e):
    u = unitroots(8)
    p1 = mpc(0.12, 0.5)
    p6 = mpc(a, 0.5)
    p7 = mpc(1.93, 0.5)
    p2 = p1 + expj(b)
    p4 = p1 + expj(c)
    p3 = p6 + expj(d)
    p5 = p3 + expj(e)
    d1 = abs(p2 - p4*u[-3]) - 1
    d2 = abs(p2 - conj(p6)*u[-1]) - 1
    d3 = abs(p3 - conj(p7)*u[2]) - 1
    d4 = abs(p4 - p7*u[1]) - 1
    d5 = abs(p5 - p5*u[1]) - 1
    vertices = symmetrise((p1, p2, p3, p4, p5, p6, p7), "D8")
    return (vertices, (d1, d2, d3, d4, d5))

@fixparams_unitdist()
def f112c():
    z1 = star_radius(28, 13)
    z2 = star_radius(28, 11) * root(1, 56, 3)
    z3 = cu(z2, z1)
    z4 = cu(0, z3, star_radius(28, 5), 1)
    vertices = symmetrise((z1, z2, z3, z4), "C28")
    return vertices

@fixparams_unitdist(-2, 1, 1.5)
def f114a(a, b, c):
    u = unitroots(19)
    p1 = mpc(a, 0.5)
    p2 = mpc(b, 0.5)
    p3 = mpc(c, 0.5)
    d1 = abs(p1 - p2*u[7]) - 1
    d2 = abs(p2 - p3*u[-1]) - 1
    d3 = abs(p3 - conj(p1)*u[-8]) - 1
    vertices = symmetrise((p1, p2, p3), "D19")
    return (vertices, (d1, d2, d3))

@fixparams_unitdist(0.72, -0.22)
def f120a(a, b):
    u = unitroots(10)
    p1 = rect(0.5, 1.646)
    p2 = -conj(p1)
    p4 = mpc(0.2425, a)
    p7 = -0.0616 + 1.128j
    p3 = cu(p1, u[2]*p7)
    p5 = cu(p1, u[1]*p4)
    p6 = cu(p2, u[-3]*p4)
    p9 = cu(p7, p2)
    p12 = cu(p6, u[-2]*p9)
    p10 = cu(0, p12, star_radius(10), 1)
    p8 = p4 + expj(b)
    p11 = cu(p8, u[-2]*p8)
    vertices = symmetrise((p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12), "C10")
    return (vertices, (abs(p5 - u[1]*p3) - 1, abs(p7 - u[2]*p11) - 1))

@fixparams_unitdist()
def f120b():
    p1 = star_radius(30, 13)
    p2 = p1 + expj(2.6895)
    p3 = cu(root(1,5,2)*p2, p2)
    p4 = cu(0, p3, star_radius(30, 7), 1)
    return symmetrise((p1, p2, p3, p4), "C30")

# The two graphs below are the remaining symmetric graph expansions
@fixparams_unitdist()
def f204a():
    p1 = 0.5289 + 1.4675j
    p2 = 1 + p1
    p3 = cu(0, p1, star_radius(34, 11), 1)
    p4 = cu(0, p1, star_radius(34, 7), 1)
    p5 = cu(0, p2, star_radius(34, 5), 1)
    p6 = cu(0, p2, star_radius(34, 3), 1)
    return symmetrise((p1, p2, p3, p4, p5, p6), "C34")

@fixparams_unitdist()
def f224c():
    p1 = 1.4425
    p2 = cu(0, p1, star_radius(56, 11), 1)
    p3 = cu(p1, 0, 1, star_radius(56, 13))
    p4 = cu(p1, 0, 1, star_radius(56, 5))
    return symmetrise((p1, p2, p3, p4), "C56")
