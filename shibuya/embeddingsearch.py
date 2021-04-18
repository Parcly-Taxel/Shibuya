from mpmath import *
from os.path import expanduser
from subprocess import run
from tempfile import NamedTemporaryFile
from functools import reduce
from shibuya.rigidity import jacobian

def symmetry_data(sym, k):
    """Given a symmetry type and index, return the GAP code defining
    that group as a finitely presented group, a mapping from generators
    to transformations of the Euclidean plane and a string describing
    allowed orbit sizes."""
    symtype, n = sym[0], int(sym[1:])
    shifts = unitroots(n, True)
    if symtype == "C":
        gap_string = f"F := FreeGroup(1);\nH := F / [F.1^{n}];"
        ops = {1: lambda z: z*shifts[k], -1: lambda z: z/shifts[k]}
        orbit_string = f"{n}"
    else:
        gap_string = f"F := FreeGroup(2);\nH := F / [F.1^{n}, F.2^2, (F.1*F.2)^2];"
        ops = {1: lambda z: z*shifts[k], -1: lambda z: z/shifts[k],
               2: conj, -2: conj}
        orbit_string = f"{n}, {2*n}"
    return (gap_string, ops, orbit_string)

def embedding_functions(edges, sym, k, gap_path):
    """Given an edge list representing a graph, a desired symmetry and
    an index describing how that symmetry gets realised in the Euclidean plane,
    return a list of embedding objects. This function depends on
    a GAP instance at gap_path and the Digraphs package there.
    sym uses Schoenflies notation, e.g. C7 or D7.
    
    An embedding object is a tuple (automorphism table, function realising embedding,
    number of variables, number of constraints)."""
    digraph_edges = [[a+1,b+1][::s] for (a, b) in edges for s in (1, -1)]
    freegroup_def, ops, orbit_sizes = symmetry_data(sym, k)
    program = f"""LoadPackage("digraphs");
Gr := DigraphByEdges({digraph_edges});
G := AutomorphismGroup(Gr);
{freegroup_def}
sH := Size(H);
for hom in IsomorphicSubgroups(G, H) do
    rings := Orbits(Image(hom));
    if not ForAll(rings, ring -> Length(ring) in [ {orbit_sizes} ]) then
        continue;
    fi;
    seeds := List(rings, ring -> Minimum(ring));
    for aut in H do
        LR := LetterRepAssocWord(Factorization(H, aut));
        block := OnTuples(seeds, aut^hom);
        Print(LR, "/", block, "\\n");
    od;
    Print("\\n");
od;
QUIT_GAP();
"""
    with NamedTemporaryFile("w+", delete=False) as f:
        f.write(program)
    proc = run([expanduser(gap_path), "-q", f.name], capture_output=True, text=True)
    res = []
    for chunk in proc.stdout.split("\n\n"):
        if not chunk:
            continue
        table = {}
        for coset in chunk.split("\n"):
            aut, points = map(eval, coset.split("/"))
            points = tuple(map(lambda v: v-1, points))
            table[tuple(aut)] = points

        ringmap = {v: i for coset in table.values() for (i, v) in enumerate(coset)}
        constraints = []
        wyckoffs = []
        for (ri, v) in enumerate(table[tuple()]):
            neigh1 = [e[1] for e in edges if e[0] == v]
            neigh2 = [e[0] for e in edges if e[1] == v]
            for w in filter(lambda x: ringmap[x] >= ri, neigh1 + neigh2):
                constraints.append((v, w))
            # Check for invariant ("Wyckoff") positions
            for g in filter(lambda aut: aut and table[aut][ri] == v, table):
                wyckoffs.append((ri, g))

        N = max(p for coset in table.values() for p in coset) + 1
        def graph_vertices(*args, table=table, constraints=constraints, wyckoffs=wyckoffs):
            vertices = [None] * N
            seeds = [mpc(*z) for z in zip(*[iter(args)]*2)]
            for (aut, coset) in table.items():
                for (ring, i) in enumerate(coset):
                    vertices[i] = reduce(lambda z, k: ops[k](z), aut, seeds[ring])
            cfuncs = [abs(vertices[i] - vertices[j])**2 - 1 for (i, j) in constraints]
            for (ri, g) in wyckoffs:
                seed = seeds[ri]
                tseed = reduce(lambda z, k: ops[k](z), g, seed)
                cfuncs.append(abs(tseed-seed)**2)
            return (vertices, cfuncs)
        res.append((table, graph_vertices, 2*len(table[tuple()]), len(constraints) + len(wyckoffs)))
    return res

def beauty_factor(G):
    """Return the "beauty factor" of an arbitrary graph, the minimum distance
    between a vertex and a non-incident edge."""
    V, E = G[0], G[1]
    dists = []
    for (i, u) in enumerate(V):
        for (j, k) in E:
            if i == j or i == k:
                continue
            v, w = V[j], V[k]
            a, b = u-v, w-v
            proj = (a.real*b.real+a.imag*b.imag) / abs(b) # scalar projection
            if 0 <= proj <= abs(b):
                dists.append(abs(a - b * proj / abs(b)))
            else:
                dists.extend((abs(a), abs(u-w)))
    return min(dists)

def embedding_run(E, sym, k, conclass, gap_path, succount=6, normlimit=1e-12,
        beautylimit=0.01, coordrange=3, maxsteps=20):
    """Generate embeddings of the given graph with the given symmetry options;
    return alongside each embedding the corresponding automorphism table and
    solution coordinates."""
    funcs = embedding_functions(E, sym, k, gap_path)
    print(f"{len(funcs)} conjugacy class(es) for {sym}:")
    for (i, quad) in enumerate(funcs):
        print(f"[{i}] -> {quad[0]} ({quad[2]} vars, {quad[3]} cons)")
    table, f, nv, nc = funcs[conclass]
    s = 0
    zfloor = eps * max(nv, nc)
    while s < succount:
        try:
            # For maximum robustness, manually implemented Newton's method
            # with pseudoinverse calculation
            x = [coordrange*(2*rand()-1) for _ in range(nv)]
            for _ in range(maxsteps):
                U, S1, V = svd(jacobian(lambda *xx: f(*xx)[1], x))
                tol = zfloor * max(S1)
                S2 = diag([1/z if z >= tol else 0 for z in S1])
                F = -matrix(f(*x)[1])
                delta = V.T * S2 * U.T * F
                x = [a + b for (a, b) in zip(x, delta)]
                if norm(delta) <= 1e-12:
                    break
            if norm(matrix(f(*x)[1])) >= normlimit:
                print("N")
                continue
            G = (f(*x)[0], E)
            if beauty_factor(G) >= beautylimit:
                yield (G, x)
                s += 1
            else:
                print("B")
        except ValueError:
            continue
