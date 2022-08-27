"""
Ed Pegg's graphs from https://math.stackexchange.com/q/3958839.
{UnitDistance, {15, 1}} is not included because it contains another
rigid graph, {UnitDistance, {12, 2}}.
"""
from mpmath import *
from shibuya.generators import cu, star_radius, ring_edges, all_unit_distances
from shibuya.graphs.rigidity import jacobian

def ud93_vertices(t):
    A = 0
    B = 1
    p0 = expj(t)
    p1 = cu(p0, B, sqrt(3), 1)
    p2 = cu(p0, A)
    p3 = cu(p2, p1)
    # ...
    p4 = cu(p0, p1)
    p5 = cu(p1, p0)
    p6 = cu(p1, p3)
    vertices = (A, B, p0, p1, p2, p3, p4, p5, p6)
    return vertices, abs(B - p3) - sqrt(3)

def ud93():
    """1 parameter, rigid"""
    f = lambda t: ud93_vertices(t)[1]
    t0 = findroot(f, [1.5, 2])
    print(diff(f, t0)) # non-zero
    vertices = ud93_vertices(t0)[0]
    return all_unit_distances(vertices)

def ud122_vertices(t, u):
    A = 0
    B = 1
    p0 = expj(t)
    p1 = 1+expj(u)
    p2 = cu(B, p1)
    p3 = cu(p0, A)
    p4 = cu(p3, p2)
    p5 = cu(p4, p3)
    p6 = cu(p2, p4)
    p7 = cu(p5, p1)
    p8 = cu(p0, p6)
    p9 = cu(p1, p0)
    vertices = (A, B, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9)
    return vertices, abs(p7 - p9) - 1, abs(p8 - p9) - 1

def ud122():
    """2 parameters, rigid"""
    f = lambda *x: ud122_vertices(*x)[1:]
    x0 = findroot(f, (1.8, 1.34))
    print(det(jacobian(f, x0))) # non-zero
    vertices = ud122_vertices(*x0)[0]
    return all_unit_distances(vertices)

def ei17_vertices(t, u):
    A = 0
    B = 1
    p0 = 1+expj(t)
    p1 = 1+expj(u)
    p2 = cu(p0, p1)
    p3 = cu(A, p0)
    p4 = cu(A, p2)
    p5 = cu(p2, A)
    p6 = cu(p1, p3)
    p7 = cu(p0, p6)
    p8 = cu(p1, p7)
    p9 = cu(p4, p8)
    p10 = cu(p9, p5)
    p11 = cu(p10, p4)
    p12 = cu(p11, p5)
    p13 = cu(p12, p9)
    p14 = cu(p13, p11)
    vertices = (A, B, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14)
    return vertices, abs(p6 - p13) - 1, abs(p8 - p14) - 1

def ei17():
    """Exoo-Ismailescu graph 17, 2 parameters, rigid"""
    f = lambda *x: ei17_vertices(*x)[1:]
    x0 = findroot(f, (2, 1.6))
    print(det(jacobian(f, x0))) # non-zero
    vertices = ei17_vertices(*x0)[0]
    return all_unit_distances(vertices)

def ei19_vertices(v):
    A = 0
    B = 1
    p0 = mpc(0.5, (sqrt(32) + sqrt(35)) / 6)
    p1 = cu(A, p0)
    p2 = cu(p0, A)
    p3 = cu(B, p0)
    p4 = cu(p0, B)
    p5 = cu(p3, p1)
    p6 = cu(p4, p2)
    p7 = cu(p5, p6)
    p8 = cu(p5, p2)
    p9 = cu(p6, p3)
    p10 = cu(p7, p0)
    # The above edges + (8, 10), (9, 10) are rigid on their own
    p11 = cu(p3, p6)
    p12 = p8 + expj(v)
    p13 = cu(p9, p12)
    p14 = cu(p11, p12)
    p15 = cu(p13, p11)
    p16 = cu(p14, p15)
    vertices = (A, B, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16)
    return vertices, abs(p9 - p16) - 1

def ei19():
    """Exoo-Ismailescu graph 19, (2, 1) parameters, rigid"""
    f = lambda v: ei19_vertices(v)[1]
    v0 = findroot(f, -0.1)
    print(diff(f, v0)) # non-zero
    vertices = ei19_vertices(v0)[0]
    return all_unit_distances(vertices)

def ei21_vertices(t, u, v, w):
    A = 0
    B = -1j
    p0 = mpc(t, u)
    p1 = cu(A, p0)
    p2 = cu(p0, A)
    p3 = cu(B, p0)
    p4 = cu(p0, B)
    p5 = cu(p4, p1)
    p6 = cu(p2, p3)
    # 1 DOF lost when (5, 6) added, but 1 DOF remains
    q1 = cu(B, p5)
    q2 = p4 + expj(v)
    q3 = cu(q1, q2)
    q4 = cu(q2, q1)
    q5 = cu(p4, q3)
    q6 = cu(q4, q5)
    r1 = cu(p6, B)
    r2 = p3 + expj(w)
    r3 = cu(r2, r1)
    r4 = cu(r1, r2)
    r5 = cu(r3, p3)
    r6 = cu(r5, r4)
    vertices = (A, B, p0, p1, p2, p3, p4, p5, p6, q1, q2, q3, q4, q5, q6, r1, r2, r3, r4, r5, r6)
    return vertices, abs(p5 - p6) - 1, abs(p2 - q6) - 1, abs(p1 - r6) - 1

def ei21(t0=0):
    """Exoo-Ismailescu graph 21, not rigid (|t0| <= sqrt(3)/2)"""
    u_estimate = 0.822878 - 0.380202*t0**2
    v_estimate = 2.7113 - 1.08206*t0
    w_estimate = 0.4303 - 1.08206*t0
    f = lambda u: ei21_vertices(t0, u, v_estimate, w_estimate)[1]
    u0 = findroot(f, u_estimate)
    f = lambda v: ei21_vertices(t0, u0, v, w_estimate)[2]
    v0 = findroot(f, [v_estimate-0.01, v_estimate+0.01])
    f = lambda w: ei21_vertices(t0, u0, v0, w)[3]
    w0 = findroot(f, [w_estimate-0.01, w_estimate+0.01])
    vertices = ei21_vertices(t0, u0, v0, w0)[0]
    return all_unit_distances(vertices)

def hodfish_vertices(t, u, v, w, x):
    A = 0
    B = 1
    p0 = 1j * expj(t)
    p1 = expj(-u)
    p2 = 1-expj(v)
    p3 = p0 + 1
    p4 = cu(p2, p1)
    p5 = cu(p4, p3)
    p6 = cu(p0, p4)
    p7 = cu(A, p5)
    p8 = cu(p6, B)
    p9 = cu(p7, p1)
    p10 = cu(p2, p8)
    q0 = p0 + expj(w)
    q1 = cu(p4, q0)
    q2 = cu(q0, p4)
    q3 = cu(B, q2)
    q4 = cu(q3, q1)
    r0 = p3 - expj(-x)
    r1 = cu(r0, p4)
    r2 = cu(p4, r0)
    r3 = cu(r2, A)
    r4 = cu(r1, r3)
    vertices = (A, B, p0, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, q0, q1, q2, q3, q4, r0, r1, r2, r3, r4)
    return vertices, abs(p0 - p10) - 1, abs(p3 - p9) - 1, abs(p0 - q4) - 1, abs(p3 - r4) - 1

def hodfish(t0=0):
    """Hochberg–O'Donnell fish graph, not rigid (|t0| <= ~0.5)"""
    u_estimate = 0.424082 - 0.311038*t0
    v_estimate = 0.424082 + 0.311038*t0
    f = lambda u, v: hodfish_vertices(t0, u, v, 0, 0)[1:3]
    u0, v0 = findroot(f, (u_estimate, v_estimate))
    w_estimate = polyval([-0.465329, 0.516923, 0.320761, 1.27248], t0)
    f = lambda w: hodfish_vertices(t0, u0, v0, w, 0)[3]
    w0 = findroot(f, w_estimate)
    x_estimate = polyval([0.465329, 0.516923, -0.320761, 1.27248], t0)
    f = lambda x: hodfish_vertices(t0, u0, v0, w0, x)[4]
    x0 = findroot(f, x_estimate)
    vertices = hodfish_vertices(t0, u0, v0, w0, x0)[0]
    return all_unit_distances(vertices)
