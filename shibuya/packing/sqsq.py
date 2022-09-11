"""
Packings of unit squares into the smallest possible unit square.
Cf. https://erich-friedman.github.io/packing/squinsqu
"""
from mpmath import *
from shibuya.draw import drawing

"""
from sympy import *
a, b, c, s, t = params = symbols("a b c s t")
R = Matrix([[1-t*t, -2*t], [2*t, 1-t*t]])
ideal = [(R @ Matrix(p))[i] - n*(1+t*t) for (p,i,n) in
         [([s-a, b-s], 0, 3), ([c-1, 1-(s-1)], 0, 2), ([a-1, s-2], 1, 2),
          ([(s-1)-2, (s-1)-1], 1, 2), ([s-(c+1), b-1], 1, 1)]]
GB = groebner(ideal, params, order="lex")
print(factor(GB[-1]).args[1])
GB = groebner(ideal, params, order="grevlex")
print(factor(GB[-2]))
ideal = [pol.subs(a, s-2) for pol in ideal]
GB = groebner(ideal, params, order="grevlex")
print("s:", solve(GB[1], s)[0])
print("c:", apart(solve(GB[0], c)[0].subs(s, solve(GB[1], s)[0])))
print("b:", solve(GB[-1], b)[0])
"""
def p11():
    t = polyroots([5, -10, -2, 14, 12, -6, 2, 2, -1])[1]
    s = (6*t+4) / ((2-t)*t+1)
    b = s - 2.5*t - 0.5/t
    c = s-1 + 1/(t+1) + 1/(t-1)
    v = mpc(2*t,1-t*t) / (1+t*t)
    r = conj(v)*1j

    k1 = mpc(s-2,s)
    k2 = mpc(s,b)
    k3 = mpc(re(r*k1)+2,im(r*mpc(s-1,s-1))) / r
    k4 = mpc(re(r*mpc(1,s-1)),im(r*mpc(1,2))) / r
    k5 = k3+k4-k1
    return (s, ((0,1), (1,2), (1j,1+1j), (c,c+1), (mpc(1,s),s*1j), (mpc(s,s),mpc(s-1,s)),
                (k1,k1-v), (k2,k2+1j*v), (k3,k3+1j*v), (k4,k4-1j*v), (k5,k5+v)))

def draw_packing(data, outfn, scale=400):
    s, sqsides = data
    res = drawing(scale, (s, s))
    for (v1, v2) in sqsides:
        v3 = v2 + 1j*(v2-v1)
        v4 = v1 + v3 - v2
        cmds = [["M", v1.real, v1.imag, v2.real, v2.imag, v3.real, v3.imag, v4.real, v4.imag], ["Z"]]
        res.add_path(cmds, {"fill": "#6dc6fb", "stroke": "#1c92cd", "stroke-width": 0.01})
    res.add_rect(0, 0, s, s, {"fill": "none", "stroke": "#000", "stroke-width": 0.01})
    res.write(outfn)
