from mpmath import unitroots

def utility():
    """Returns the minimal integral embedding of K_{3,3} with maximum edge length 2.
    This graph is rigid, but not first-order/infinitesimally rigid."""
    vertices = unitroots(6)
    edges = ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (0, 3), (1, 4), (2, 5))
    return (vertices, edges)
