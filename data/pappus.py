# Draws a unit-distance embedding of the Pappus graph.

u12 = unitroots(12)
inner = [u / 2 for u in u12[1::2]]
outer_s = cu(inner[1], inner[5])
outer = [outer_s * u for u in u12[::2]]
tilt_s = cu(0, outer[5])
tilt = [tilt_s * u for u in u12[::2]]
vertices = inner + outer + tilt
edges = [(0, 3), (1, 4), (2, 5)] + [(i, (i + 1) % 6 + 6) for i in range(6)] + [(i, (i - 1) % 6 + 6) for i in range(6)] + \
        [(i + 12, (i + 1) % 6 + 12) for i in range(6)] + [(i + 6, (i + 1) % 6 + 12) for i in range(6)]
