#!/usr/bin/env python3
from sympy import *
a, b, c, d, s, t, x = symbols('a b c d s t x')

# The printout of this program is a Singular program which finds the
# degree-32 minimal polynomial of the t parameter of the unit-distance
# embedding of the F26A graph.

def ul(a, b):
    v = a - b
    return minpoly(v.dot(v), x).subs(x, 1)

def pcp(a, b):
    return a[0]*b[1]-a[1]*b[0]

ct = sqrt(1-t**2)
r60 = Matrix([[1, -s], [s, 1]]) / 2
p1 = Matrix([t, ct]) / 2
p1r = r60 @ p1
ph = Matrix([-1, 0])
p2 = Matrix([a, b])
E1 = ul(p2, ph)
E2 = ul(p2, p1r)
p3 = Matrix([c, d])
E3 = ul(p3, p1)
E4 = ul(p3, p2)
print('LIB "elim.lib";')
print("option(prot);")
print("ring r=(0,s),(a,b,c,d,t),dp; minpoly=s2-3;")
print(f"poly p1={factor(E1)};")
print(f"poly p2={E2.subs(s**2, 3) / -16};")
print(f"poly p3={factor(E3)};")
print(f"poly p4={factor(E4)};")
t1 = r60 @ (p1 + p2 - p3)
t2 = -p3
t3 = r60.T @ p3
u1 = t1 - t2
u2 = t3 - t2
u3 = t1 - t3
K = 4*pcp(u1, u2)**2 - u1.dot(u1) * u2.dot(u2) * u3.dot(u3)
K = K.expand().collect(ct)
Kwith, Kwithout = K.coeff(ct, 1), K.coeff(ct, 0)
K = Kwithout**2 - Kwith**2 * (1-t**2)
K = factor(K.subs(s, sqrt(3)))
print(f"poly p5={K.args[-1]};")
print("ideal I=p1,p2,p3,p4,p5;")
print("ideal J=elim(I,1..4);")
print("factorize(J[1]);")
print("quit;")
