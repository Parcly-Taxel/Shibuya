"""
Packings of points (or equivalently circles) into a unit circle
maximising the minimum separation d. Each function below returns
a pair (d, sequence of points).
Cf. https://erich-friedman.github.io/packing/cirincir
and http://hydra.nat.uni-magdeburg.de/packing/cci/cci.html
"""
from mpmath import *
from shibuya.draw import drawing
from shibuya.generators import cu

def thas(t):
    """Given t = tan(x/2), return the point on the unit circle
    (cos(x), sin(x))."""
    return mpc((1-t)*(1+t),2*t) / (1+t*t)

def psmall(n):
    if n > 9:
        return None
    centre = n > 6
    m = n-centre
    d = 2*sinpi(fraction(1,m))
    points = unitroots(m)
    if centre:
        points.append(0)
    return (d, points)

"""
option(prot); LIB "elim.lib";
ring r=(0,s),(v,u,t,a,d),dp; minpoly=s2-3;
ideal I=(1-t2)-(s*d/2+a)*(1+t2),4t-d*(1+t2),
        d2*(1+u2)-4u2,(1-v2-(a-d)*(1+v2))^2+4v2-d2*(1+v2)^2,
        v*(1-(3tu+3u2)+tu3)-((t+3u)-(u3+3tu2));
I=sat(I,d*(d2-2)*(d6+(s-10)*d4+(-7s+32)*d2+(12s-31)))[1];
poly pd=elim(I,tuva)[1];
pd;
"""
def p10():
    s = sqrt(3)
    up = [[0, 1], [-2, -2], [2, 1], [8, 8], [-8, 11], [-6, -6], [6, -5]]
    ap = [[0, 1], [1, -14], [-12, 76], [53, -201], [-104, 265], [82, -151], [-8, 14]]
    u = polyroots([polyval(z, s) for z in up], extraprec=50)[0]
    d = 2*u / hypot(1,u)
    a = sqrt(polyroots([polyval(z, s) for z in ap], extraprec=50)[1])
    zu = thas(u)
    return (d, [zu**(k/2) for k in range(-7,8,2)] + [a, a-d])

def p11():
    d = 2*sinpi(fraction(1,9))
    points = unitroots(9)
    z = mpc(re(2*points[3]-points[4]),im(points[4]))
    points += [z, conj(z)]
    return (d, points)

"""
option(prot); LIB "elim.lib";
ring r=(0,s),(u,d),dp; minpoly=s2-3;
ideal I=d2*(1+u2)-4u2,((1-u2)-(1+u2)*d/(2s))^2+(2u-(1+u2)*d/2)^2-d2*(1+u2)^2;
I=sat(I,d+s)[1];
poly pd=elim(I,u)[1];
pd;
"""
def p12():
    s = sqrt(3)
    d = polyroots([1, -5, 7, 0, -27, 27])[1] / s
    zu = thas(sqrt((2/d)**2 - 1))
    points = [u*z for u in unitroots(3) for z in (-1, d/s, zu, conj(zu))]
    return (d, points)

def p13():
    d = 2*sinpi(fraction(1,10))
    points = unitroots(10) + [0.35, -0.22+0.32j, -0.22-0.32j]
    return (d, points)

"""
option(prot); LIB "elim.lib";
ring r=(0,s),(w,v,u,t,c,b,a,d),dp; minpoly=s2-3;
ideal I=d2*(1+t2)-4t2,(1-u2)-(s*d/2+a)*(1+u2),2u-d/2*(1+u2),
        c2+((a-b)/2)^2-d2,
        ( ((a+b)/2*(1+v2)-(1-v2))^2+(c*(1+v2)-(2v))^2-(d*(1+v2))^2 ) / (1+v2),
        ( (b*(1+w2)-(1-w2))^2+(2w)^2-(d*(1+w2))^2 ) / (1+w2),
        v*(1-(t2+2tu))-((2t+u)-t2u),
        w*(1-(t2+2tv))-((2t+v)-t2v);
poly pd=elim(I,wvutcba)[1];
"""
def p14():
    tp = [1240029, 8148762, 8444007, -47475396, -98474049, 107193618, 346191165, -99412272, -613361646, 3630420, 586703574, 27276264, -310427154, -4021164, 90370890, -1870128, -13879647, 233874, 1143027, 5436, -48189, -870, 841]
    with extraprec(100):
        s = sqrt(3)
        t = polyroots(tp, extraprec=100)[4] * s
        d = 2*t / hypot(1,t)
        zt = thas(t)
        a = re(sqrt(zt)) - s*d/2
        c = cu(a, zt**2.5, d, d)
        b = 2*re(c) - a
    return (d, [zt**(k/2) for k in range(-9,10,2)] + [a, b, c, conj(c)])

def p15():
    h = 0.5 / tan(pi/5) + 1
    d = 1 / hypot(h,0.5)
    z = mpc(h,0.5) * d
    points = [u*y for u in unitroots(5) for y in (z, conj(z), z-d)]
    return (d, points)

def p16():
    pass # TODO

def p17():
    pass # TODO

def p18():
    d = 2*sinpi(fraction(1,12))
    points = unitroots(12)
    z = 1 + 1j*(1-points[-1])
    points += [z*u for u in unitroots(6)]
    return (d, points)

def p19():
    dat = p18()
    return (dat[0], dat[1] + [0])

def pn(n):
    """Return the best known point packing in a unit circle
    with the given number of points."""
    if n < 10:
        return psmall(n)
    return eval(f"p{n}()")

def draw_packing(data, outfn, scale=400):
    d, points = data
    res = drawing(scale, (2+d, 2+d), (-(1+d/2), -(1+d/2)))
    res.add_circle(0, 0, 1, {"fill": "none", "stroke": "#000", "stroke-width": 0.005*d})
    for p in points:
        res.add_circle(p.real, p.imag, d/2, {"fill": "#6dc6fb", "fill-opacity": "0.8",
                                             "stroke": "#1c92cd", "stroke-width": 0.005*d})
        res.add_circle(p.real, p.imag, 0.02*d)
    res.add_circle(0, 0, 1+d/2, {"fill": "none", "stroke": "#000", "stroke-width": 0.005*d})
    res.write(outfn)
