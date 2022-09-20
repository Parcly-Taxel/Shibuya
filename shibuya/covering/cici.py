"""
The canonical disc covering problem: minimise r in a covering of the
unit circle by N circles of radius r. Each function below returns
a pair (r, sequence of centres).
Cf. https://erich-friedman.github.io/packing/circovcir
"""
from mpmath import *
from shibuya.draw import drawing
from shibuya.generators import cu

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
        t = cu(-1, a-r, r, r)
        m = (cu(0, a, 1, r) + cu(t, 0, r, 1)) / 2
    return (r, [a, t, m, conj(t), conj(m)])

"""
LIB "elim.lib";
ring r=0,(E,D,C,B,e,d,c,b,a,q),dp;
ideal I=b^2+c^2-q,B^2+C^2-q,d^2+e^2-q,D^2+E^2-q,
(a+2b+d)^2+e^2-1,(a+B)^2+C^2-1,(a+B+2D)^2+(C+2E)^2-1,
(B+2D-b-d)^2+(C+2E-c-e)^2-q;
poly p=factorize(elim(I,BCDEedc)[1])[1][4];
poly p1=factorize(resultant(p,diff(p,a),b))[1][3];
poly p2=factorize(resultant(p,diff(p,b),b))[1][3];
poly q_p=factorize(resultant(p1,p2,a))[1][5];
q_p;
poly a_p=factorize(resultant(p2,q_p,q))[1][3];
a_p;
"""
def c6():
    with extraprec(100):
        qp = [7841367, -33449976, 62607492, -63156942, 41451480, -19376280, 5156603, -746832, 54016, 3072]
        ap = [7841367, -90811125, 344734812, -376618374, -189418782, 373328980, 71968935, -130579181, 27056736, -800000]
        bp = [7841367, -7610085, -71679924, -45430974, 58085154, 26617188, -12289505, -5546269, 2814464, -327680]
        r = sqrt(polyroots(qp)[1])
        a = -sqrt(polyroots(ap)[1])
        b = sqrt(polyroots(bp)[0])
        t = cu(a, b, r, r) + cu(0, b, 1, r) - b
        m = (cu(a, 0, r, 1) + cu(0, t, 1, r)) / 2
    return (r, [a, b, t, m, conj(t), conj(m)])

def c_rad(n):
    m = n-1
    r = 1 / (1 + 2*cospi(fraction(2,m)))
    d = 2*cospi(fraction(1,m)) * r
    return (r, [0] + [d*u for u in unitroots(m)])

def c11():
    with extraprec(100):
        rp = [16, -56, 121, -458, 1709, -1773, -4406, 12770, -4789, -16653, 19746, -5805, -12708, 15406, -5897, 729]
        ap = [704, 1040, 13239, 47050, 572272, 940099, 3697162, 42769344, 27697007, 178794966, 268623532, 23638520, -77077232, -1705056, 3613248, -279936]
        r = polyroots(rp)[2]
        a = polyroots(ap)[3]
        p1 = cu(a+r, 1, r, r)
        t = cu(0, p1, 1, r)
        p2, v = t-r, t-p1
        p3 = cu(0, p2, 1, r) + a+r - p1
        b = 2*re(p3 - v) - a
        c = -re(v) - sqrt(1 - im(v)**2)
        p4 = cu(c, b, r, r) + v
        ps = [p1, p2, p3, p4]
    return (r, [a, b, c] + ps + [conj(p) for p in ps])

def c12():
    o = unitroots(3)
    r = polyroots([1, -1, 3, -1])[0]
    z = 1 + r*o[1]
    return (r, [p*u for p in (1-2*r, 2-1/r, z, conj(z)) for u in o])

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
