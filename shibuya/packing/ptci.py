"""
Packings of points (or equivalently circles) into a unit circle
maximising the minimum separation d. Each function below returns
a pair (d, sequence of points).
Cf. https://erich-friedman.github.io/packing/cirincir
and http://hydra.nat.uni-magdeburg.de/packing/cci/cci.html
"""
from mpmath import *
from shibuya.draw import drawing

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

def p14():
    pass # TODO

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
