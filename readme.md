**Shibuya** is a utility for drawing graphs, with a particular focus on integral embeddings in the plane (which every graph has) and unit-distance embeddings (which not all graphs have).

This project started in my [Dounreay](https://gitlab.com/parclytaxel/Dounreay) "incubator" repository, initially just for exploration into integral graphs. I managed to find a unit-distance embedding of the Foster graph, but that was essentially all I found out on my own. After a question by Ed Pegg regarding the _rigidity_ of a certain unit-distance graph – [I answered it, showing it was rigid](https://math.stackexchange.com/a/3955860/357390) – I realised that keeping the relevant files in Dounreay would be too much of a hassle, so I split them off here.

The name recalls that of another of my repositories dealing heavily with networks of some kind, [Shinjuku](https://gitlab.com/parclytaxel/Shinjuku). Shibuya Station is the second-busiest train station in the world, behind Shinjuku.

# The rendering system

How Shibuya renders graphs is very similar to yet another project of mine, [Malibu](https://gitlab.com/parclytaxel/Malibu); it also depends on [mpmath](http://mpmath.org).

The key function is `draw(graph, outfn, scale, vertsize, edgewidth)`:

* `graph` is a pair of `(vertices, edges)`, where `vertices` is a sequence of complex numbers and `edges` is a sequence of pairs of indices into `vertices`.
* `outfn` is the output file name (`.svg` is automatically appended).
* One unit in raw graph coordinates corresponds to `scale` pixels in the SVG file.
* `vertsize` and `edgewidth` are the radius of each vertex and the width of each edge in raw graph coordinates respectively.

The last three arguments are optional.
