import xml.etree.ElementTree as ET

class drawing:
    def __init__(self, scale=400, canvas=None, offset=None):
        self.scale = scale
        ET.register_namespace("", "http://www.w3.org/2000/svg")
        if canvas is None: # (width, height)
            canvas = (1, 1)
        if offset is None: # (x, y)
            offset = (0, 0)
        self.root = ET.Element("svg", {"width": str(canvas[0]*scale), "height": str(canvas[1]*scale),
                "viewBox": f"{offset[0]*scale} {-(offset[1]+canvas[1])*scale} {canvas[0]*scale} {canvas[1]*scale}"})
        self.tree = ET.ElementTree(self.root)

    def add_object(self, tag, mains, styledict=None):
        attrib = mains.copy()
        if styledict is not None:
            if "stroke-width" in styledict:
                styledict["stroke-width"] *= self.scale
            attrib["style"] = ";".join(f"{k}:{v}" for (k,v) in styledict.items())
        self.root.append(ET.Element(tag, attrib))

    def add_circle(self, x, y, r=0.02, styledict=None):
        scale = self.scale
        attrib = {"cx": str(x*scale), "cy": str(-y*scale), "r": str(r*scale)}
        self.add_object("circle", attrib, styledict)

    def add_path(self, cmds, styledict=None):
        """Accepts a sequence of sequences of a letter followed by numbers
        representing SVG path commands."""
        scale = self.scale
        cmdstrs = []
        for cmd in cmds:
            head, rest = cmd[0], cmd[1:]
            rest = [(-1)**i * x*scale for (i, x) in enumerate(rest)]
            cmdstrs.append(head + " ".join(map(str, rest)))
        attrib = {"d": "".join(cmdstrs)}
        self.add_object("path", attrib, styledict)

    def add_edge(self, x1, y1, x2, y2, styledict=None):
        if styledict is None:
            styledict = {"fill": "none", "stroke": "#000", "stroke-width": 0.005}
        self.add_path([("M", x1, y1, x2, y2)], styledict)

    def write(self, fn):
        self.tree.write(f"{fn}.svg", "unicode")

def draw_graph(graph, outfn, scale=400, pad=0.04):
    """Draw graph, specified as (vertices, edges), writing to outfn.svg.
    One unit in graph's coordinates corresponds to scale pixels in the output;
    a padding of pad units is applied all around."""
    vertices, edges = graph
    reals = [v.real for v in vertices]
    imags = [v.imag for v in vertices]
    x = min(reals) - pad
    y = min(imags) - pad
    width = max(reals) - min(reals) + 2*pad
    height = max(imags) - min(imags) + 2*pad
    res = drawing(scale, (width, height), (x, y))

    for (i1, i2) in edges:
        v1 = vertices[i1]
        v2 = vertices[i2]
        res.add_edge(v1.real, v1.imag, v2.real, v2.imag)
    for v in vertices:
        res.add_circle(v.real, v.imag)
    res.write(outfn)

def draw_pointpacking(data, outfn, scale=400):
    """Draw the given packing of points in a unit square, writing to outfn.svg."""
    d, points = data
    res = drawing(scale)
    for p in points:
        res.add_circle(p.real, p.imag, d/2, {"fill": "#6dc6fb", "fill-opacity": "0.8",
                                             "stroke": "#1c92cd", "stroke-width": 0.005*d})
        res.add_circle(p.real, p.imag, 0.02*d)
    res.write(outfn)
