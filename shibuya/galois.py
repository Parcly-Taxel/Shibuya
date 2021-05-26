"""
Tools for finding Galois groups of polynomials with coefficients in algebraic extensions,
with a view to getting exact expressions for the coordinates that appear in
Shibuya's embeddings.
"""
from mpmath import *
from itertools import combinations, permutations

# Orbit-length partitions from Soicher & McKay '85 (https://core.ac.uk/download/pdf/81182361.pdf)
# with respect to the modes in resolvent(). The leftmost column gives the LMFDB index;
# "+" in the d column means the quadratic resolvent polynomial is reducible.
"""
deg 3:
  d 2-
1 + 3,3
2   6

deg 4:
  d 2+  2-
1   2,4 4^3
2 + 2^3 4^3
3   2,4 4,8
4 + 6   12
5   6   12

deg 5:
  d 2+  2-   e5
1 + 5^2 5^4
2 + 5^2 10^2
3   10  20   5,10
4 + 10  20
5   10  20   15

deg 6:
   d 2+    3+      2-     e6
1    3,6^2 2,6^3   6^5
2    3^3,6 2,6^3   6^5
3    3,6^2 2,6,12  6,12^2
4  + 3,12  4^2,6^2 6,12^2
5    6,9   2,18    6^2,18
6    3,12  6^2,8   6,12^2
7  + 3,12  4^2,12  6,24
8    3,12  8,12    6,24   8^2,12^2
9    6,9   2,18    12,18  2^2,18^2
10 + 6,9   2,18    12,18
11   3,12  8,12    6,24   8^2,24
12 + 15    10^2    30
13   6,9   2,18    12,18  4,36
14   15    20      30     20^2
15 + 15    20      30
16   15    20      30     40

deg 7:
  d 2+  3+     2-
1 + 7^3 7^5    7^6
2   7^3 7^3,14 14^3
3 + 21  7^2,21 21^2
4   21  14,21  42
5 + 21  7,28   42
6 + 21  35     42
7   21  35     42
"""

def resolvent(mpt, i, in_pol, mode="d", prec=500,
        digitlimit=100, maxsteps=1000, choptol=30):
    """Print Singular code to factorise the resolvent polynomial of in_pol,
    a univariate polynomial with coefficients in QQ[t], given as a coefficient list
    of polynomial coefficient lists (the inner lists are evaluated at t = polyroots(mpt)[i]).
    The resolvent invariant to use is controlled by mode;
    prec, digitlimit, maxsteps and choptol control the working precision."""
    with extradps(prec):
        extdeg = len(mpt) - 1 # degree of the extension QQ(t)
        mpt_str = "".join(f"{c:+}t{extdeg-i}" for (i, c) in enumerate(mpt) if c != 0)
        mpt_str = mpt_str.lstrip("+")
        t = polyroots(mpt)[i]
        basis = [t**n for n in range(extdeg)]

        inroots = polyroots([polyval(coeff, t) for coeff in in_pol], extraprec=mp.dps)
        if mode == "2+":
            resroots = [a+b for (a,b) in combinations(inroots, 2)]
        elif mode == "3+":
            resroots = [a+b+c for (a,b,c) in combinations(inroots, 3)]
        elif mode == "2-":
            resroots = [a-b for (a,b) in permutations(inroots, 2)]
        elif mode == "d":
            sqrtdisc = fprod(a-b for (a,b) in combinations(inroots, 2))
            resroots = [-sqrtdisc, sqrtdisc]
        elif mode == "e5":
            hadamard = ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2))
            resroots = [(comb[i1] + comb[i2] - comb[i3] - comb[i4])**2
                    for comb in combinations(inroots, 4) for (i1,i2,i3,i4) in hadamard]
        elif mode == "e6":
            sqrtdisc = fprod(a-b for (a,b) in combinations(inroots, 2))
            resroots = [a+b+c+sign*sqrtdisc for (a,b,c) in combinations(inroots, 3) for sign in (-1, 1)]
        y_deg = len(resroots)
        respol = [1] + resroots
        for d in range(1, y_deg+1):
            r = respol[d]
            respol[d] = 0
            for i in range(d, 0, -1):
                respol[i] -= r * respol[i-1]
        resmonos = [f"y{y_deg}"]
        for (r, c) in enumerate(respol[1:], 1):
            if (Z := chop(c, tol=mpf(f"1e-{mp.dps-choptol}"))) == 0:
                continue
            bZ = basis + [Z]
            coords = pslq(bZ, maxcoeff=10**digitlimit, maxsteps=maxsteps)
            print(r, coords, nstr(fdot(bZ, coords)))
            if coords[-1] > 0:
                coeff_str = "".join(f"{-c:+}/{coords[-1]}t{i}" for (i, c) in enumerate(coords[:-1]))
            else:
                coeff_str = "".join(f"{c:+}/{-coords[-1]}t{i}" for (i, c) in enumerate(coords[:-1]))
            coeff_str = coeff_str.lstrip("+")
            resmonos.append(f"+({coeff_str})*y{y_deg-r}")
        res_str = "".join(resmonos)
        print(f"ring r=(0,t),y,dp; minpoly={mpt_str};")
        print(f"poly p={res_str};")
        print("factorize(p);")
