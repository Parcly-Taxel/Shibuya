# Draws a unit-distance graph containing a regular heptagon that was conjectured by
# Ed Pegg to be rigid in https://math.stackexchange.com/q/3954719/357390
# I proved it to be rigid, even when vertices 7 and 14 in the realisation below are deleted

s1 = subtend_r(1, 7)
r1 = [s1 * u for u in unitroots(7)]
s2 = subtend_r(2, 7)
s2 = cc(0, s1, s2, 1)
r2 = [s2 * u for u in unitroots(7)]
s3 = subtend_r(3, 7)
s3 = cc(s1, 0, 1, s3)
r3 = [s3 * u for u in unitroots(7)]
vertices = r1 + r2 + r3
edges = ring_edges(7, ((0, 0, 1), (1, 1, 2), (0, 1, 0), (2, 2, 3), (0, 2, 0), (1, 2, 0)))

A = zeros(len(edges), 2*len(vertices))
for (r, (i1, i2)) in enumerate(edges):
    v1, v2 = vertices[i1], vertices[i2]
    dreal, dimag = re(v1 - v2), im(v1 - v2)
    A[r,2*i1] = dreal
    A[r,2*i1+1] = dimag
    A[r,2*i2] = -dreal
    A[r,2*i2+1] = -dimag
