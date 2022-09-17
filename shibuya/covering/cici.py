"""
The canonical disc covering problem: minimise r in a covering of the
unit circle by N circles of radius r. Each function below returns
a pair (r, sequence of centres).
Cf. https://erich-friedman.github.io/packing/circovcir
"""
from mpmath import *
from shibuya.draw import drawing

def thassqrt(q):
    return mpc(1-q,2*sqrt(q)) / (1+q)

def c3():
    return (sqrt(3)/2, [0.5*u for u in unitroots(3)])

def c4():
    return (sqrt(0.5), [sqrt(0.5)*u for u in unitroots(4)])

"""
option(prot); LIB "elim.lib";
ring r=0,(t,u,c,b,a,s),dp;
ideal I=2b+1+s-a,
        (a*(1+t2) - (1-t2))^2 + (2t)^2 - s2*(1+t2)^2,
        (b*(1+u2) - (1-u2))^2 + (c*(1+u2) - (2u))^2 - s2*(1+u2)^2,
        (b+1)^2 + c2 - s2,
        ((1-t2)*(1+u2) - (1-u2)*(1+t2))^2 + (2t*(1+u2) - 2u*(1+t2))^2 - 4s2*(1+t2)^2*(1+u2)^2;
I=I[1],I[2]/(1+t2),I[3]/(1+u2),I[4],I[5]/((1+t2)*(1+u2));
poly p=elim(I,tubc)[1]/s;
ideal J=p,diff(p,a);
factorize(elim(J,a)[1]);
"""
def c5():
    with extraprec(50):
        r = polyroots([1296, 2112, -3480, 1360, 1665, -1776, 22, -800, 625])[2]
        a = polyroots([1296, -4576, -2232, 2552, 4261, -568, 7664, -3328, -1024])[1]
        z1 = mpc(a-r-1, sqrt((r+a+1) * (3*r-a-1))) / 2 # mpc(b,c)
        t2 = polyroots([100, -3304, 39337, -158014, 287575, -223724, 33895, -622, -75])[1]
        u2 = polyroots([675, 918, -207, -5300, -5759, -250, 1007, 824, -100])[3]
        z2 = (thassqrt(t2)+thassqrt(u2)) / 2
    return (r, [a, z1, conj(z1), z2, conj(z2)])

def c_rad(n):
    m = n-1
    r = 1 / (1 + 2*cospi(fraction(2,m)))
    d = 2*cospi(fraction(1,m)) * r
    return (r, [0] + [d*u for u in unitroots(m)])

def cn(n):
    """Return the best known circle covering of a unit circle
    with the given number of circles."""
    if 7 <= n <= 10:
        return c_rad(n)
    return eval(f"c{n}()")

def draw_packing(data, outfn, scale=400):
    r, centres = data
    res = drawing(scale, (2,2), (-1,-1))
    res.add_circle(0, 0, 1.008)
    res.add_circle(0, 0, 1, {"fill": "#fff"})
    for c in centres:
        res.add_circle(c.real, c.imag, r, {"opacity": "0.1"})
    res.write(outfn)
