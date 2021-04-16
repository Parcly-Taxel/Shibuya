from mpmath import re, im, zeros, diff

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
