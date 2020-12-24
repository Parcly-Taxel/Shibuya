from mpmath import unitroots
from shibuya.generators import star_radius, ring_edges

def circulant(n, taps):
    """Return the circulant graph on n vertices with offsets given by taps.
    The vertices are placed on an n-gon of unit side length."""
    r = star_radius(n)
    vertices = [r*u for u in unitroots(n)]
    edges = ring_edges(n, [(0, 0, k) for k in taps])
    return (vertices, edges)

def cycle(n):
    """Return the cycle graph on n vertices."""
    return circulant(n, [1])

def complete(n):
    """Return the complete graph on n vertices."""
    return circulant(n, range(1, n//2+1))

def mobiusladder(n=3):
    """Return the Möbius ladder on 2n vertices. The n = 3 case corresponds
    to the utility graph; its embedding returned by this function is integral
    with edge lengths 1 and 2, as well as rigid, but not first-order rigid."""
    return circulant(2*n, (1, n))
