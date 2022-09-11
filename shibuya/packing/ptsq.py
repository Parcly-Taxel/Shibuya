"""
Packings of points (or equivalently circles) into a unit square
maximising the minimum separation d. Each function below returns
a pair (d, sequence of points).
Cf. https://erich-friedman.github.io/packing/cirinsqu
and http://hydra.nat.uni-magdeburg.de/packing/csq/csq.html
"""
from mpmath import *
from shibuya.draw import drawing

def chickenwire(a, b):
    """Return the chicken-wire point packing where the points are placed
    at alternate vertices of an a-by-b-rectangle grid, stretched to fill
    the unit square."""
    return (sqrt(fraction(a*a+b*b, a*a*b*b)),
            [mpc(fraction(x,a), fraction(y,b)) for x in range(a+1) for y in range(x%2,b+1,2)])

def mappacking(a, b):
    """Return the map-like point packing where an a-by-b point lattice is
    folded like a map to fit in the unit square. Circles contact on the
    b (vertical) axis."""
    a_, b_ = a-1, b-1
    d = (a_*b_ - sqrt(a_*a_ - b_*b_ + 1))/(a_*(b_*b_ - 1))
    y = 1-b_*d
    return (d, [mpc(fraction(i,a_), i%2*y+j*d) for i in range(a) for j in range(b)])

def p2():
    return (sqrt(2), [0, 1+1j])

def p3():
    return (sqrt(6)-sqrt(2), [0, mpc(1,2-sqrt(3)), mpc(2-sqrt(3),1)])

def p4():
    return (1, [0, 1, 1j, 1+1j])

def p5():
    return (sqrt(0.5), [0, 1, 1j, 1+1j, 0.5+0.5j])

def p6():
    return chickenwire(2, 3)

def p7():
    d = 4-2*sqrt(3)
    return (d, [0, d, d*1j, d*(1+1j), mpc(d/2,1), mpc(1,d/2), 0.95+0.95j])

def p8():
    return (sqrt(2 - sqrt(3)), [bp*z + 0.5+0.5j for z in unitroots(4)
                                                for bp in (0.5+0.5j, (1-sqrt(3))/2)])

def p9():
    return (0.5, [mpc(x,y)/2 for x in (0,1,2) for y in (0,1,2)])

"""
from sympy import *
params = symbols("a b y x d")
ideal = (S("x^2+y^2-d^2"),
         S("(2*x+d-1)^2+b^2-d^2"),
         S("a^2+(1-2*y)^2-d^2"),
         S("(a+d-1)^2+(b+d-1)^2-d^2"),
         S("(x-a-d)^2+(y-1+d)^2-d^2"))
GB = groebner(ideal, params, order="lex")
print(GB)
"""
def p10():
    dp = [1180129, -11436428, 98015844, -462103584, 1145811528,
          -1398966480, 227573920, 1526909568, -1038261808, -2960321792,
          7803109440, -9722063488, 7918461504, -4564076288, 1899131648,
          -563649536, 114038784, -14172160, 819200]
    xp = [4833808384, -18292893696, 231089016320, -1126322897408, 2641176117632,
          -3604472915840, 3142380993216, -1893047603328, 926983988992, -469985292736,
          225222946304, -56978206656, -15097884104, 17327716152, -5236535740,
          336548136, 165929282, -38942194, 2584905]
    d = polyroots(dp, extraprec=50)[0]
    x = polyroots(xp, extraprec=50)[4]
    y = sqrt(d*d - x*x)
    a = sqrt(d*d - (1-2*y)**2)
    b = sqrt((1-2*x) * (2*(d+x)-1))
    return (d, [mpc(*z) for z in ((0,0), (x,y), (2*x,0), (0,2*y), (2*x+d,0),
                                  (a,1), (a+d,1), (1,b), (1,b+d), (a+d,1-d))])

def p11():
    d = polyroots([1, 8, -22, 20, 18, -24, -24, 32, -8])[1]
    x = polyroots([16, -32, 16, 8, -7])[1] * d
    y = polyroots([16, -32, 16, -8, 1])[1] * d
    return (d, [mpc(*z) for z in ((0,0), (0,d), (d,0), (1,1), (1-sqrt(2)*d,1), (1,1-sqrt(2)*d),
                                  (1-d/sqrt(2),1-d/sqrt(2)), (x,y), (y,x), (0.91,0.025), (0.025,0.91))])

def p12():
    return chickenwire(3, 5)

"""
option(prot); LIB "elim.lib";
ring r=(0,s),(b,a,y,x,d),dp; minpoly=s2-3;
ideal I=x2+y2-d2,a2+b2-d2,(2a+d-1)^2+(2y+d-1)^2-d2,(2x-1)^2+(2b+d-1)^2-d2,
        (x-(1-(a+s*b)/2))^2+(y-(1-2*b+(b-s*a)/2))^2-d2;
I=sat(I,d)[1];
poly pd=elim(I,abxy)[1];
pd;
"""
def p13():
    dp = [[0, 864665116641711], [-28104584929111908, 47659330746396264], [-645768516385207920, 1114660079628188624],
          [-271002221507689328, 475367780231855840], [149768958170972320, -252188349821109840], [16368378056297492608, -28371625172339499520],
          [-40884498213224701312, 70824275651796366336], [26436355874572015104, -45777262685375092736], [43275912192017374208, -74972437684294657792],
          [-113150947454306898944, 195987208674775394304], [123266390798296981504, -213498628394498887680], [-76447923840235352064, 132407451118567415808],
          [17244750471993139200, -29868073120995438592], [21010669645192298496, -36390914668133941248], [-31976204717970653184, 55384018697947971584],
          [25006610638992310272, -43312662998256189440], [-13114829814314827776, 22715578710630858752], [4706963012829839360, -8152714665959358464],
          [-1109734524459155456, 1922120126119804928], [155129707824676864, -268692944792322048], [-9771872514211840, 16925399418142720]]
    xp = [[0, 73784756620092672], [539842596297395328, -1605913295417536896], [7649548279968554880, -10234806945572200896],
          [-101918184555224364000, 167525102239286754144], [641333830739022614528, -1090702875060607724832], [-2659372602272398688672, 4570572090402518188512],
          [7771488504603038564480, -13409472436692001524240], [-16663159117033256164168, 28801139619747100254312], [27107560866385811343448, -46893072115504231885992],
          [-34325683395193177037508, 59406798153660548536188], [34404368061616303059628, -59558925568511328332472], [-27509287558718541027184, 47630405389431326784288],
          [17527456045223510456866, -30350741803402788860364], [-8798107223686951373501, 15235922140930525178403], [3391589658274816259101, -5873551105602317787456],
          [-956451470050696868585, 1656417047045976187689], [180345709940951615817, -312329185811717220222], [-19464222788927823607, 33707487455411514249],
          [1395775840746124655, -2416989309515807532], [-360023509218944039, 623542079396343255], [61976449083692325, -107345194733251650]]
    ap = [[0, 2235901715760384], [60574205995728768, -126082435118172288], [-1696161125425460352, 3030870107487473472],
          [11698123941573183456, -20515150730375641376], [-45340880818061884512, 79018331582022385280], [82216046888637149408, -143109156853286167968],
          [-72880964682792091824, 127056231762824542112], [22972231413740749480, -40582382208463716040], [12919760502136984680, -21729310160910330360],
          [-14629683999098483188, 24883952247436017500], [2925984654080788552, -4791097672773392692], [4102943269153386896, -7252409486230078200],
          [-3589076397294650268, 6282906018324080674], [617723127797077571, -1095894459412061093], [492621788137825195, -844629291514062038],
          [-159505178683494885, 273883093024406604], [-60459706781338571, 105257244763333305], [32046567838768767, -55600559431419026],
          [-6644850582641145, 11521181021601304], [2830038050144638, -4902736797205390], [-637851253724384, 1104827687429172]]
    with extraprec(50):
        s = sqrt(3)
        d = polyroots([polyval(z, s) for z in dp], extraprec=50)[0]
        x = polyroots([polyval(z, s) for z in xp], extraprec=50)[0]
        a = polyroots([polyval(z, s) for z in ap], extraprec=50)[0]
        y = sqrt(d*d - x*x)
        b = sqrt(d*d - a*a)
        q = (1-(a+s*b)/2, 1-2*b+(b-s*a)/2)
        return (d, [mpc(*z) for z in ((0,0), (x,y), (2*x,0), (0,2*y), (0,2*y+d), q, (0.37,0.64),
                                      (1,1), (1-a,1-b), (1-2*a,1), (1,1-2*b), (1-2*a-d,1), (1,1-2*b-d))])

def p14():
    d = 2*(4-sqrt(3))/13
    return (d, [mpc(x,y)*d for x in range(3) for y in range(3)] +
               [mpc(*z) for z in ((d*0.5,1), (1,d*0.5), (d*1.5,1), (1,d*1.5), (0.965,0.965))])

def p15():
    d = (1+sqrt(2)-sqrt(3))/2
    cd75 = 1 - d*(sqrt(6)+sqrt(2))/4
    d15 = d*(sqrt(6)-sqrt(2))/4
    points = []
    for x in range(4):
        for y in range(4):
            if x == y == 2:
                continue
            if x == 2:
                px = cd75
                py = d15+d*y if y < 3 else 1-d15
            elif y == 2:
                py = cd75
                px = d15+d*x if x < 3 else 1-d15
            else:
                px = x*d if x < 3 else 1
                py = y*d if y < 3 else 1
            points.append(mpc(px,py))
    return (d, points)

def p16():
    return (mpf(1)/3, [mpc(x,y)/3 for x in range(4) for y in range(4)])

"""
option(prot); LIB "elim.lib";
ring r=0,(a,b,k,d),dp;
ideal I=(d-1/2)^2+a^2-d^2,(a-1/2)^2+(d+1/2-b)^2-d^2,
        (k-1/2)^2+(a/2-b/2)^2-d^2,k^2+(a/2+b/2-d)^2-d^2;
I=sat(I,d*(d-2))[1];
poly pd=elim(I,abk)[1];
pd;
"""
def p17():
    d = polyroots([1, -4, 6, -14, 22, -20, 36, -26, 5])[0]
    a = sqrt(d-0.25)
    b = d+0.5 - sqrt((d-a+0.5) * (d+a-0.5))
    k = sqrt((a+b) * (4*d-a-b)) / 2
    points = [mpc(*z) for z in ((0,0), (d,0), (0,d), (0,2*d), (0,1), (a,d+0.5), (k,(a+b)/2))]
    points = points + [1-conj(z) for z in points] + [mpc(0.5,a), mpc(0.5,b), mpc(0.5,0.99)]
    return (d, points)

def p18():
    return chickenwire(4, 6)

"""
option(prot); LIB "elim.lib";
ring r=(0,s),(a,y,x,d),dp; minpoly=s2-3;
ideal I=x2+y2-d2,(x-a)^2+(d+3y-1)^2-d2,
        (1-(s/2+1)*d-a)^2+(3d/2+2y-1)^2-d2,2x+s*d-1;
I=sat(I,d*(2d+1-s))[1];
poly pd=elim(I,ayx)[1];
pd;
"""
def p19():
    s = sqrt(3)
    dp = [[0, 22], [25, -65], [-280, -422], [213, 467], [-78, -162], [52, -52]]
    d = polyroots([polyval(z, s) for z in dp], extraprec=50)[1]
    x = (1-s*d)/2
    y = sqrt(d*d - x*x)
    dx, dy = d*s/2, d/2
    a = x + sqrt((1-3*y) * (2*d+3*y-1))
    return (d, [0, 1] + [mpc(*z) for z in ((0,0), (dx,dy), (2*dx,0), (1-x,y), (1,0),
                                           (0,d), (x,d+y), (0.5,y+dy), (1-dx,2*y+dy), (1,2*y),
                                           (0,d+2*y), (x,d+3*y), (2*x,d+2*y), (1-dx,1.5*d+2*y), (1,d+2*y),
                                           (a,1), (a+d,1), (0.03,0.96), (0.97,0.96))])

def p20():
    return mappacking(5, 4)

def p25():
    return (0.25, [mpc(x,y)/4 for x in range(5) for y in range(5)])

def p27():
    return chickenwire(5, 8)

def p30():
    return mappacking(6, 5)

def pn(n):
    """Return the best known point packing in a unit square
    with the given number of points."""
    return eval(f"p{n}()")

def draw_packing(data, outfn, scale=400):
    d, points = data
    res = drawing(scale, (1+d, 1+d), (-d/2, -d/2))
    res.add_rect(0, 0, 1, 1, {"fill": "none", "stroke": "#000", "stroke-width": 0.005*d})
    for p in points:
        res.add_circle(p.real, p.imag, d/2, {"fill": "#6dc6fb", "fill-opacity": "0.8",
                                             "stroke": "#1c92cd", "stroke-width": 0.005*d})
        res.add_circle(p.real, p.imag, 0.02*d)
    res.add_rect(-d/2, -d/2, 1+d, 1+d, {"fill": "none", "stroke": "#000", "stroke-width": 0.005*d})
    res.write(outfn)
