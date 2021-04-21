import numpy as np
from math import gcd
from os.path import expanduser
from subprocess import run
from tempfile import NamedTemporaryFile
rng = np.random.default_rng()

conclasses_template = """LoadPackage("digraphs");
Gr := DigraphByEdges({0});
G := AutomorphismGroup(Gr);
H := {1};

for hom in IsomorphicSubgroups(G, H) do
    Ih := Image(hom);
    rings := Orbits(Ih);
    if not ForAll(rings, ring -> Length(ring) in [ {2} ]) or
            Number(rings, ring -> Length(ring) = 1) > 1 then
        continue;
    fi;
    seeds := List(rings, Minimum);

    aut_table := [ ];
    for aut in H do
        word := ExtRepOfObj(Factorization(H, aut));
        block := OnTuples(seeds, aut^hom);
        Add(aut_table, [ word, block ]);
    od;

    GetWords := function(elms)
        local preims;
        preims := List(elms, g -> PreImagesRepresentative(hom, g));
        return List(preims, g -> ExtRepOfObj(Factorization(H, g)));
    end;

    ring_descs := [ ];
    for i in [1..Length(seeds)] do
        seed := seeds[i];
        # Get distinct neighbours up to H's action
        stab := Stabilizer(Ih, seed);
        neigh_orbits := Orbits(stab, OutNeighborsOfVertex(Gr, seed));
        distinct_neighs := List(neigh_orbits, Minimum);
        neigh_rings := List(distinct_neighs, w -> PositionProperty(rings, ring -> w in ring));
        forward_indices := PositionsProperty(neigh_rings, x -> i <= x);
        # For each forward neighbour find the H-action sending the corresponding seed to it
        sources := seeds{{neigh_rings{{forward_indices}}}};
        destinations := distinct_neighs{{forward_indices}};
        permuters := ListN(sources, destinations, {{x,y}} -> RepresentativeAction(Ih, x, y));
        neigh_words := GetWords(permuters);
        # Get transformations leaving the seed invariant
        invariant_words := GetWords(GeneratorsOfGroup(stab));
        ring_desc := [ invariant_words, neigh_rings{{forward_indices}}, neigh_words ];
        Add(ring_descs, ring_desc);
    od;

    Print([ aut_table, ring_descs ], "\\n\\n");
od;
QUIT_GAP();
"""

def symmetry_data(sym):
    """Given a symmetry type, return the GAP code defining that group,
    allowed orbit sizes (orders of site symmetry groups of Wyckoff positions)
    and a function mapping indexes to realisations of the group in Euclidean space
    (functions mapping external representations to orthogonal NumPy matrices)."""
    symtype, n = sym[0], int(sym[1:])
    cyclic_shifts = list(filter(lambda m: gcd(n,m) == 1, range(1, n//2+1)))
    if symtype == "C":
        groupdef = f"CyclicGroup(IsPermGroup, {n})"
        orbit_sizes = [1, n]
    else:
        groupdef = f"DihedralGroup(IsPermGroup, {2*n})"
        orbit_sizes = [1, n, 2*n]
    def rindexer(k):
        def rfunc(extrep):
            A = np.eye(2)
            base_theta = 2*cyclic_shifts[k]*np.pi/n
            for (g, times) in zip(*[iter(extrep)] * 2):
                if g == 1:
                    s, c = np.sin(times*base_theta), np.cos(times*base_theta)
                    A = np.array([[c, -s], [s, c]]) @ A
                else:
                    A = np.diag([1, -1]) @ A
            return A
        return rfunc
    return (groupdef, orbit_sizes, rindexer)

def embedding_conclasses(edges, sym, gap_path):
    """Given a graph's edge list and a desired symmetry, return conjugacy classes
    of graph embeddings respecting said symmetry. This function depends on
    a GAP instance at gap_path and the Digraphs package there.
    sym uses Schoenflies notation, e.g. C7 or D7."""
    digraph_edges = [[a+1,b+1][::s] for (a, b) in edges for s in (1, -1)]
    groupdef, orbit_sizes, rfunc = symmetry_data(sym)
    orbit_sizes = ", ".join(map(str, orbit_sizes))
    program = conclasses_template.format(digraph_edges, groupdef, orbit_sizes)
    with NamedTemporaryFile("w+", delete=False) as f:
        f.write(program)
    proc = run([expanduser(gap_path), "-q", f.name], capture_output=True, text=True)
    chunks = []
    for chunk in proc.stdout.split("\n\n"):
        if not chunk:
            continue
        chunk = eval(chunk)
        chunk.append(rfunc)
        chunks.append(chunk)
    return chunks

def conclass_realisation(chunk, k):
    """Realise the given conjugacy class (chunk) according to the given index."""
    aut_table, ring_descs, rindexer = chunk
    rfunc = rindexer(k)
    I = rfunc([])
    N = I.shape[0]
    coord_mats = [np.zeros((N,0))]
    for (invariant_words, _, _) in ring_descs:
        if not invariant_words:
            coord_mats.append(np.eye(N))
        else:
            A = np.concatenate([rfunc(w) - np.eye(N) for w in invariant_words])
            _, Sigma, VT = np.linalg.svd(A)
            tol = N*len(invariant_words) * np.finfo(float).eps * max(Sigma)
            coord_mats.append(VT[sum(Sigma > tol):].T)
    starts = np.cumsum([C.shape[1] for C in coord_mats])
    nv = starts[-1]
    strides = [C.shape[1] for C in coord_mats]
    constraints = []
    for (ri, (xx, neigh_rings, neigh_words)) in enumerate(ring_descs, 1):
        for (rj, word) in zip(neigh_rings, neigh_words):
            C = np.zeros((N, nv))
            C[:,starts[ri-1]:starts[ri-1]+strides[ri]] += coord_mats[ri]
            C[:,starts[rj-1]:starts[rj-1]+strides[rj]] -= rfunc(word) @ coord_mats[rj]
            constraints.append(C)
    preF = np.stack(constraints) # @ this with x and take squared norm minus one to get F(x)
    Tmats = np.stack([rfunc(pair[0]) for pair in aut_table]) # transformation matrices
    Tverts = np.stack([pair[1] for pair in aut_table]) # vertices corresponding to each Tmat
    return (preF, Tmats, Tverts, np.concatenate(coord_mats, axis=1), starts)

def save_embeddings(E, symtype, gap_path):
    for (cc, chunk) in enumerate(embedding_conclasses(E, symtype, gap_path)):
        k = 0
        while 1:
            try:
                preF, Tmats, Tverts, coord_mat, starts = conclass_realisation(chunk, k)
                if k == 0:
                    print(f"{symtype}, cc = {cc}: nv = {preF.shape[2]}, nc = {preF.shape[0]}")
                np.savez(f"{symtype}-{cc}-{k}", preF, Tmats, Tverts, coord_mat, starts)
                k += 1
            except IndexError:
                break

def load_test_embedding(symtype, cc, k, successes=np.inf, failures=np.inf,
        coordrange=3, maxsteps=20):
    with np.load(f"{symtype}-{cc}-{k}.npz") as arrd:
        preF, Tmats, Tverts, coord_mat, starts = arrd.values()
        nv = preF.shape[2]
        strides = [starts[i+1] - starts[i] for i in range(len(starts)-1)]
    preJ = 2 * np.swapaxes(preF, 1, 2) @ preF

    s, f = 0, 0
    while s < successes and f < failures:
        x = rng.uniform(-coordrange, coordrange, nv)
        for _ in range(maxsteps):
            c = preF @ x
            F = (c*c).sum(axis=1) - 1
            J = preJ @ x
            delta = np.linalg.lstsq(J, -F, rcond=None)[0]
            x += delta
            if np.linalg.norm(delta) <= 1e-12:
                break
        c = preF @ x
        F = (c*c).sum(axis=1) - 1
        if np.linalg.norm(F) > 1e-12:
            f += 1
            continue
        s += 1
        vertices = [None] * Tverts.max()
        for (i, M) in enumerate(Tmats):
            for (j, v) in enumerate(Tverts[i]):
                sl = slice(starts[j], starts[j]+strides[j])
                coords = x[sl]
                C = coord_mat[:,sl]
                if vertices[v-1] is None:
                    vertices[v-1] = M @ C @ coords @ [1, 1j]
        yield (x, vertices)

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
