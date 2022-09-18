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
        z1xp = [335994734583, -925055302392, 1245681478080, -1737278445780, 1649211336432, -1213918522608, 517124784064, -136608531264, 17765003264, -663552000]
        z1yp = [335994734583, -665599287231, 1897433047512, -2672527493250, 2452227909282, -430846892340, -903975402593, 814316672985, -232138895616, 7464960000]
        z4xp = [282289212, -664669152, 1485276912, -2162565189, 2039760072, -1482840408, 758989690, -225747200, 20917924, -253125]
        z4yp = [31365468, -117592236, 439727724, -870853095, 1390690842, 1186401779, -230544808, -235749633, 12250566, 11390625]
        r = sqrt(polyroots(qp)[1])
        a = -sqrt(polyroots(ap)[1])
        b = sqrt(polyroots(bp)[0])
        z1x = -sqrt(polyroots(z1xp)[1])
        z1y = sqrt(polyroots(z1yp)[2])
        z1 = mpc(z1x,z1y)
        z4x = sqrt(polyroots(z4xp)[1])
        z4y = sqrt(polyroots(z4yp)[2])
        z4 = mpc(z4x,z4y)
    return (r, [a, b, z1, z4, conj(z1), conj(z4)])

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
