"""
Cubic symmetric graphs. Most of the embeddings realised here were taken from MathWorld.
"""
from mpmath import *
from shibuya.generators import cu, star_radius, ring_edges, all_unit_distances

# F4A = tetrahedron() or complete(4) (not unit-distance)
# F6A = circulant(6, (1, 3)) or mobiusladder(3) (not unit-distance)
# F8A = genpetersen("cube")
# F10A = genpetersen("petersen")

def heawood():
    """Return the symmetric unit-distance embedding of the Heawood graph (F14A)
    tucked away in Mathematica's GraphData."""
    P = [10485760, 78643200, 263192576, 543686656, 812777472, 942080000, 843317248, 552468480, 208879616, -31170560, -99213312, -76779520, -32795648, 7878144, 17269760, 16256512, 11392032, 4836080, 3014064, 361320, 69498, -165789]
    # v0 is the only real root of the above polynomial
    v0 = polyroots(P, maxsteps=1000)[0]
    p0 = mpc(0.5, v0)
    p1 = mpc(sqrt(1-(v0+0.5)**2)-0.5, -0.5)
    p2 = cu(p0, -p0)
    p3 = cu(p1, -p1)
    p4 = cu(p2, p3)
    vertices = [mpc(s*re(v), im(v)) for s in (1, -1) for v in (p0, -p0, p1, -p1, p2, p3, p4)]
    return all_unit_distances(vertices)

# F16A = genpetersen("mobiuskantor")

def pappus():
    """Return a unit-distance embedding of the Pappus graph (F18A)."""
    u6 = unitroots(6)
    r0 = [u*0.5j for u in u6]
    z1 = cu(r0[2], r0[0])
    r1 = [z1*u for u in u6]
    z2 = cu(0, z1)
    r2 = [z2*u for u in u6]
    vertices = r0 + r1 + r2
    edges = ring_edges(6, ((0, 0, 3), (0, 1, 0), (0, 1, -2), (2, 2, 1), (1, 2, 0)))
    return (vertices, edges)

# F20A = genpetersen("dodecahedron")
# F20B = genpetersen("desargues")
# F24A = genpetersen("nauru")

def f26a_vertices(t):
    A, B, C = unitroots(6)[4:1:-1]
    p2 = mpc(t, sqrt(1-t**2)) / 2
    p1 = p2 * root(1, 6, 1)
    p3 = p2 * root(1, 6, 5)
    p4 = cu(p1, B)
    p5 = cu(p2, C)
    p6 = cu(p3, -A)
    p7 = cu(p1, -p6)
    p8 = cu(p4, p2) # 1
    p9 = cu(p5, p3) # 1
    p10 = cu(-p8, p7) # 1
    V = (A, B, C, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
    return [v*s for v in V for s in (1, -1)]

def f26a():
    """Return a unit-distance embedding of the F26A graph."""
    L = [225763882172416, 677291646517248, -168428484689920, 20022957145325568, 110425351460487168, -62271608227627008, -1121399379324829696, 109067800210833408, 7540122093258539008, -7773427560423096320, -68122575204409933824, 52458094831308111872, 538023325894893109248, 78711684450812166144, -2822235439482864205824, -2995029967245104644096, 8042320072275722354688, 17883673052584642527232, -4091384557677847359488, -44565059159935463940096, -34337101102504632257536, 36701690560501134167040, 70580403666521815358208, 14232955138192862263296, -47035732072645104214528, -35808690608913191615616, 5894307940195308883552, 16637521343183451078624, 4349659847836042273980, -2275416792867611843748, -1127656410109393840797, 66655446601201742988, 72730432421368293597]
    with extradps(10):
        t = findroot(lambda x: polyval(L, x), 0.4) / 2
    return all_unit_distances(f26a_vertices(t))

def coxeter():
    """Return a unit-distance embedding of the Coxeter graph (F28A)."""
    u7 = unitroots(7)
    s1 = star_radius(7)
    s2 = star_radius(7, 2)
    s3 = star_radius(7, 3)
    r0 = [-s2*u for u in u7]
    r1 = [s3*u for u in u7]
    z2 = cu(r0[0], r1[3])
    r2 = [z2*u for u in u7]
    z3 = cu(0, z2, s1, 1)
    r3 = [z3*u for u in u7]
    vertices = r0 + r1 + r2 + r3
    edges = ring_edges(7, ((0, 0, 2), (1, 1, 3), (3, 3, 1), (0, 2, 0), (1, 2, -3), (2, 3, 0)))
    return (vertices, edges)

def tutte8_vertices(x):
    u5 = unitroots(5)
    r0 = [u*0.25j for u in u5]
    r1 = [-u*x*1j for u in u5]
    z2 = cu(r0[2], r0[4])
    r2 = [z2*u for u in u5]
    z3 = cu(r0[3], r1[1])
    r3 = [z3*u for u in u5]
    z4 = cu(r1[4], r1[3])
    r4 = [z4*u for u in u5]
    z5 = cu(r3[0], r2[0])
    r5 = [z5*u for u in u5]
    return (r0 + r1 + r2 + r3 + r4 + r5, abs(z4-z5)-1)

def tutte8():
    """Return a unit-distance embedding of the Tutte 8-cage (F30A; MathWorld calls this
    the Levi graph)."""
    x0 = findroot(lambda x: tutte8_vertices(x)[1], [0.33, 0.34])
    vertices = tutte8_vertices(x0)[0]
    edges = ring_edges(5, ((0, 2, 1), (0, 2, 3), (0, 3, 2), (1, 3, 4), (1, 4, 1), (1, 4, 2),
                           (2, 5, 0), (3, 5, 0), (4, 5, 0)))
    return (vertices, edges)

def dyck():
    """Return a unit-distance embedding of the Dyck graph (F32A)."""
    r0 = unitroots(8)
    r1 = [sqrt(2)*u for u in r0]
    z2 = cu(r0[1], 0, 1, star_radius(8))
    r2 = [z2*u for u in r0]
    z3 = cu(0, r1[0], star_radius(8, 3), 1)
    r3 = [z3*u for u in r0]
    vertices = r0 + r1 + r2 + r3
    edges = ring_edges(8, ((0, 1, 1), (0, 1, -1), (0, 2, -1), (2, 2, 1), (1, 3, 0), (3, 3, 3)))
    return (vertices, edges)

def klein(a1=4.47, a2=2.42, a3=0.7, s1=1, s2=-1):
    """Return a unit-distance embedding of the cubic Klein graph (F056B)."""
    u7 = unitroots(7)
    z0 = star_radius(7)
    r0 = [z0*u for u in u7]
    z1 = z0 + expj(a1)
    z2 = z1 + expj(a2)
    z3 = z1 + expj(a3)
    r1 = [z1*u for u in u7]
    r2 = [z2*u for u in u7]
    r3 = [z3*u for u in u7]
    z4 = cu(*(r2[2], r3[0])[::s1])
    z5 = cu(*(r2[0], r3[1])[::s2])
    r4 = [z4*u for u in u7]
    r5 = [z5*u for u in u7]
    z6 = cu(0, r4[0], star_radius(7, 2), 1)
    z7 = cu(0, r5[0], star_radius(7, 3), 1)
    r6 = [z6*u for u in u7]
    r7 = [z7*u for u in u7]
    vertices = r0 + r1 + r2 + r3 + r4 + r5 + r6 + r7
    edges = ring_edges(7, ((0, 0, 1), (0, 1, 0), (1, 2, 0), (1, 3, 0),
                           (2, 4, -2), (3, 4, 0), (2, 5, 0), (3, 5, -1),
                           (4, 6, 0), (5, 7, 0), (6, 6, 2), (7, 7, 3)))
    return (vertices, edges)

def foster_vertices(r):
    v3a = 0.265
    v3 = v3a * root(1, 5, 2)
    v2 = cu(v3, v3a)
    v5r = root(1, 20, 7) * r
    v5r2 = -v5r.conjugate()
    v5 = v5r * root(1, 15, 14)
    v0 = cu(v5r, v5r2)
    v1 = cu(v2, v0)
    v4 = cu(v3, v5)
    vgens = (v0, v1, v2, v3, v4, v5)
    vertices = [v*u for v in vgens for u in unitroots(15)]
    return (vertices, abs(v1 - v4*root(1, 15, 2)) - 1)

def foster():
    """Return a unit-distance embedding of the Foster graph (F90A)."""
    r0 = findroot(lambda r: foster_vertices(r)[1], 0.35)
    vertices = foster_vertices(r0)[0]
    edges = ring_edges(15, ((0, 1, 0), (1, 2, 0), (2, 3, 0), (3, 4, 0), (4, 5, 0), (5, 0, -1),
                            (0, 5, -2), (2, 3, -6), (4, 1, -2)))
    return (vertices, edges)

def biggssmith():
    """Return a unit-distance embedding of the Biggs–Smith graph (F102A)."""
    s1 = star_radius(17)
    s2 = star_radius(17, 2)
    s4 = star_radius(17, 4)
    s8 = star_radius(17, 8)
    u17 = unitroots(17)
    r1 = [s1*u*1j for u in u17]
    r4 = [s4*u*1j for u in u17]
    r8 = [s8*u*-1j for u in u17]
    sh1 = cu(r1[0], r4[0])
    rh1 = [sh1*u for u in u17]
    sh2 = cu(sh1, r8[7])
    rh2 = [sh2*u for u in u17]
    s2 = cu(sh2, 0, 1, s2)
    r2 = [s2*u for u in u17]
    vertices = r1 + r4 + rh1 + r8 + rh2 + r2
    edges = ring_edges(17, ((0, 0, 1), (1, 1, 4), (3, 3, 8), (5, 5, 2),
                            (0, 2, 0), (1, 2, 0), (2, 4, 0), (4, 5, 0), (4, 3, 7)))
    return (vertices, edges)

def heawood_gerbracht():
    """Return a unit-distance embedding of the Heawood graph.
    This is the first construction given in Gerbracht (2008),
    Eleven Unit Distance Embeddings of the Heawood Graph,
    https://arxiv.org/abs/0912.5395"""
    p5, l5, p7, l7, p2, l3 = 0, 1, 1+1j, 1+2j, 2j, 1j
    # Degree-79 polynomial for the x-coordinate of l4
    P = [82521703002365615643033600000, 152135800369825007098920960000, -2120259444356145889512456192000, -8175821639408563679884718899200, 11025799477301561380923592949760, 149189048927171391219263917572096, 341989984727973884867396338188288, -763800345871643605733535512788992, -6892489761595983453459595854256128, -14303114368662112977785429692643328, 31343179682405215504161837658819584, 254616663098419271111012531383618560, 477056183905245971917488031692938304, -898635822877066299154282762314323520, -6556557400356413683063078157405200320, -11463391775661584618715895715025904128, 14201705397119143149709337683063717104, 109385892925207478360122518287948266224, 200727376265061817580032667984094835280, -60617631026953339799378305296984824656, -1097279690260575519531876572540059803892, -2435243231716218466580115477132980137292, -1650827959998751884275421145646879272940, 4733784662326174469816987234959768253776, 17321709733106215547946139735151891780269, 28955348159492426037443729536713509636773, 21867253654523569285667250014704999794577, -30934416501269415569285918492882277029311, -152272756904971138353148344210050406803617, -325157218431048323421805399697113970403121, -451645674349824891937650532097542435080453, -343872926425618669220741202688368202345065, 286935408276107233753158822122577885606822, 1923952833473734147443634652898764867278198, 5180867575272248126071836905848828341927070, 9840137643451726574992603743314811193317886, 12720991312674958659494390034544200285598942, 6140751881298455069763046326781936407849238, -19570606574427556470966233073766236628787234, -67869160289243415287139367956058055810404822, -126213399593210124294769323126921742027688497, -164022007275895644197052849670790737036540873, -146045321267662575006252144965793560225509061, -57662664820854923809493690824194000968797973, 77130998985650655864689962089382720577858101, 213238754173051016042819729417269617854966165, 355269696471069385886754716351566237266830009, 664103288660372783854070699409333594554864741, 1533983381070251025995839971747580678500964852, 3550298683130683434462662943215234037891507412, 7274584518541872070335933586256322019748139172, 12902691890291653798206974719870993995735753540, 19954407479150801176386566760213350973570196080, 27201188778043412156622512508710379868716241416, 32963773017875955980864755706102737727974961688, 35716781564909427260214641236767872783162088204, 34722813139066795200081139797717329223025992699, 30346538554876120431728853077314517314609386819, 23863989284324858347511498529857889181323950967, 16888459659695355326863471817880692818622047623, 10751995223766688842173817330681019107518783545, 6152912915312070842952691100801887803370907305, 3160933625571584072448347845721693351127774301, 1455265549140319863369871581645012065857955441, 599083193386406195758633497190777431543886358, 219897164806211674807756610736580167553542758, 71715126275516155874072784490774971066237326, 20690863770430719393270631202371992513434414, 5253121604626527413065008160498494678879110, 1166012291532956694933924468283307736346382, 224472408717775611491021156521892619843158, 37109973679879574898320679050599920287450, 5203227805425306398124203036880713293545, 608930205226991194133708856923335926849, 58239553681851019741523172701651095197, 4422420653730204080254904433581059629, 255652807673380729611728470237761555, 10528063279784456967398200502468691, 273675328487397647237991825000783, 3348011046054687446588586894387]
    with extradps(100):
        xl4 = findroot(lambda x: polyval(P, x), [-0.74, -0.72], maxsteps=100)
    yl4 = sqrt(4 - (1-xl4)**2)
    l4 = mpc(xl4, yl4)
    p4 = (l4 + l5) / 2
    p6 = cu(l4, l7)
    p3 = cu(l3, l4)
    l2 = cu(p4, p2)
    l1 = p3 + 1
    l6 = cu(p5, p6)
    p1 = cu(l6, l1)
    vertices = (l1, l2, l3, l4, l5, l6, l7, p1, p2, p3, p4, p5, p6, p7)
    edges = ring_edges(7, ((0, 1, 0), (0, 1, -1), (0, 1, 2)))
    return (vertices, edges)
