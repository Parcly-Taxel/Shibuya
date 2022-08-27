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

def hessian(f, x0):
    """Construct the Hessian matrix of the R^n -> R function f at x0."""
    n = len(x0)
    H = zeros(n,n)
    for i in range(n):
        for j in range(n):
            dvec = tuple(int(k == i) + int(k == j) for k in range(n))
            H[i,j] = diff(f, x0, dvec)
    return H

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
