#!/usr/bin/env python3
import sys
from mpmath import *
mp.dps = 100

def cu(z1, z2):
    """Constructs the point at a distance of 1 (hence Construct Unit) from z1 and z2,
    left of the line from z1 to z2. Both arguments are mpc's."""
    m = (z1 + z2) / 2
    d, theta = polar(z2 - z1)
    return m + sqrt(1 - d * d / 4) * expj(theta + pi / 2)

def cc(z1, z2, r1, r2):
    """Constructs the point at a distance of r1 from z1 and r2 from z2, left of the
    line from z1 to z2 (hence Construct Cosine (rule))."""
    d, theta0 = polar(z2 - z1)
    theta = acos((d * d + r1 * r1 - r2 * r2) / (2 * d * r1))
    return z1 + r1 * expj(theta0 + theta)

def subtend_r(p, q):
    """Calculates the radius of the circle if a unit-length chord on it subtends
    p/q of a full rotation. This is used to precisely place points on a star of edges."""
    return sqrt(0.5 / (1 - cospi(2 * mpf(p) / mpf(q))))

@extradps(50)
def newton(pol_list, x0):
    """Uses Newton's method to find a root of the polynomial given by pol_list,
    where the highest-degree term is first. Starts from x0.
    This is used when mpmath's polyroots fails to find the desired root."""
    delta, x = 1, mpf(x0)
    while abs(delta) > mpf(f"1e-{mp.dps - 30}"):
        p_at, d_at = polyval(pol_list, x, True)
        delta = p_at / d_at
        x -= delta
    return x

def ring_edges(N, triples):
    """Suppose the vertices are grouped into some number of rings, each ring having N vertices.
    Those vertices can then be grouped into N congruent units.
    This function returns an index list for the edges constructed from triples,
    where a triple (a, b, k) means "a vertex in ring a is connected to the vertex in ring b, k units forward."."""
    L = []
    for (a, b, k) in triples:
        L.extend([(i + a * N, (i + k) % N + b * N) for i in range(N)])
    return L

# Define default parameters
scale = 400
uscale = 1
# Eval the given file, helped by the functions above, to obtain vertices and edges.
# .py is automatically added onto the end.
gname = sys.argv[1]
with open(f"data/{gname}.py", 'r') as f:
    exec(f.read())

def draw_graph(vertices, edges, scale):
    # Convert vertices to ordinary complex numbers
    vertices = [complex(chop(v) / uscale) for v in vertices]
    # Calculate drawing parameters
    vr = 0.02 * scale
    esw = 0.005 * scale
    reals = [v.real for v in vertices]
    imags = [v.imag for v in vertices]
    x = (min(reals) - 0.05) * scale
    y = (-max(imags) - 0.05) * scale
    width = (max(reals) - min(reals) + 0.1) * scale
    height = (max(imags) - min(imags) + 0.1) * scale
    # Now scale all vertices
    vertices = [v * scale for v in vertices]
    
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="{x} {y} {width} {height}">']
    pathedges = []
    for (i1, i2) in edges:
        v1 = vertices[i1]
        v2 = vertices[i2]
        pathedges.append(f'M{v1.real} {-v1.imag} {v2.real} {-v2.imag}')
    pathedges = "".join(pathedges)
    lines.append(f'  <path d="{pathedges}" stroke="#000" stroke-width="{esw}"/>')
    
    for v in vertices:
        lines.append(f'  <circle cx="{v.real}" cy="{-v.imag}" r="{vr}"/>')
    lines.append('</svg>')
    
    with open(f"{gname}.svg", 'w') as f:
        print('\n'.join(lines), file=f)

draw_graph(vertices, edges, scale)
