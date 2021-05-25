from mpmath import *

def rigidity_matrix(graph):
    """Return the rigidity matrix of the given graph.
    If this matrix has rank 2*|V|-3 the graph is infinitesimally rigid."""
    vertices, edges = graph
    A = zeros(len(edges), 2*len(vertices))
    for (r, (i1, i2)) in enumerate(edges):
        v1, v2 = vertices[i1], vertices[i2]
        dreal, dimag = re(v1 - v2), im(v1 - v2)
        A[r,2*i1] = dreal
        A[r,2*i1+1] = dimag
        A[r,2*i2] = -dreal
        A[r,2*i2+1] = -dimag
    return A

def jacobian(f, x0):
    """Construct the Jacobian matrix of the (possibly multivariate)
    function f at x0."""
    n, m = len(x0), len(f(*x0))
    J = zeros(m,n)
    for i in range(m):
        fi = lambda *x: f(*x)[i]
        for j in range(n):
            dvec = tuple(int(k == j) for k in range(n))
            J[i,j] = diff(fi, x0, dvec)
    return J

def findroot_svd(f, x0, maxsteps=20, normlimit=1e-12):
    """Find a root of f using Newton's method starting from x0,
    using the SVD to solve the linear system for robustness.
    maxsteps has the same meaning as in findroot(); the final
    result must satisfy norm(f(x*))^2 <= normlimit to be accepted."""
    x = list(x0)
    n, m = len(x), len(f(*x))
    zfloor = eps * max(m,n)
    for _ in range(maxsteps):
        U, S1, V = svd(jacobian(f, x))
        tol = zfloor * max(S1)
        S2 = diag([1/z if z >= tol else 0 for z in S1])
        F = -matrix(f(*x))
        delta = V.T * S2 * U.T * F
        x = [a + b for (a, b) in zip(x, delta)]
        if norm(delta) <= 1e-12:
            break
    fx = f(*x)
    if fdot(fx, fx) >= normlimit:
        raise ValueError(f"no convergence in {maxsteps} steps")
    return x

def resolvent():
    pass

"""from itertools import combinations, permutations
mp.dps = 500
from sympy import symbols, primitive
tt, y = symbols('t y')

t = tan(pi/5)
roots1 = polyroots([200000, -340000*t**3 + 3300000*t - 800000, 1241000*t**3 - 1350000*t**2 - 12075000*t + 13650000, -7348000*t**3 + 4499000*t**2 + 69876000*t - 42035000, 17785750*t**3 - 14172000*t**2 - 167986250*t + 131909000, -35837175*t**3 + 25526125*t**2 + 338456625*t - 239087875, 48037075*t**3 - 35278500*t**2 - 454618625*t + 333445000, -47826050*t**3 + 34426200*t**2 + 453559450*t - 328578000, 34254150*t**3 - 24618850*t**2 - 325370850*t + 236356350, -17669250*t**3 + 12736700*t**2 + 167836650*t - 121882300, 6468960*t**3 - 4701465*t**2 - 61045600*t + 44354575, -1647255*t**3 + 1228155*t**2 + 15239325*t - 11088125, 278315*t**3 - 229280*t**2 - 2487865*t + 1821840, -28860*t**3 + 24200*t**2 + 237140*t - 174000, 1324*t**3 - 1264*t**2 - 10020*t + 7440], maxsteps=1000, extraprec=500)
sqrtdisc = fprod(a-b for (a,b) in permutations(roots1, 2))
resultant_roots = [sqrtdisc, -sqrtdisc]
# resultant_roots = [a+b for (a,b) in combinations(roots1, 2)]
ydeg = len(resultant_roots)
P = y**ydeg
for r in range(1, ydeg+1):
    Z = chop(fsum(fprod(x) for x in combinations(resultant_roots, r)), tol=1e-300)
    if Z == 0:
        continue
    basis = [1, t, t**2, t**3, Z]
    coords = pslq(basis, maxcoeff=10**50, maxsteps=1000)
    print(r, coords, nstr(fdot(basis, coords)))
    Z_expr = sum(c*tt**i / -coords[-1] for (i, c) in enumerate(coords[:-1]))
    P += (-1)**r * Z_expr * y**(ydeg-r)
P = primitive(P)[1]
print("ring r=(0,t),y,dp; minpoly=t4-10t2+5;")
print(f"poly p={P};")
print("factorize(p);")"""
