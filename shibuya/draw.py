from mpmath import chop

def draw(graph, outfn, scale=400, vertsize=0.02, edgewidth=0.005, vertstyle=None):
    """Draw the given graph (specified as (vertices, edges)),
    writing to the given filename with .svg automatically appended.
    One unit in the raw graph coordinates corresponds to scale pixels in the output;
    the vertices are vertsize units in radius and the edges are edgewidth units in width.
    
    If vertstyle is provided, it must be a function mapping vertex numbers to SVG style strings;
    returning None or the empty string keeps the vertex black."""
    vertices, edges = graph
    verts = [complex(chop(v) * scale) for v in vertices]
    reals = [v.real for v in verts]
    imags = [v.imag for v in verts]
    pad = 2 * vertsize * scale
    x = min(reals) - pad
    y = -max(imags) - pad
    width = max(reals) - min(reals) + 2 * pad
    height = max(imags) - min(imags) + 2 * pad
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="{x} {y} {width} {height}">']

    pathdata = []
    for (i1, i2) in edges:
        v1 = verts[i1]
        v2 = verts[i2]
        pathdata.append(f'M{v1.real} {-v1.imag} {v2.real} {-v2.imag}')
    lines.append(f'  <path d="{"".join(pathdata)}" stroke="#000" stroke-width="{edgewidth * scale}"/>')
    for (vi, v) in enumerate(verts):
        style = ""
        if vertstyle and (props := vertstyle(vi)):
            style = f' style="{props}"'
        lines.append(f'  <circle cx="{v.real}" cy="{-v.imag}" r="{vertsize * scale}"{style}/>')
    lines.append('</svg>')

    with open(f"{outfn}.svg", 'w') as f:
        print(*lines, sep='\n', file=f)
