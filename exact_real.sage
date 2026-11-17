"""Exact check that the 14 points P_k = (1, 1 + t_k^2, t_k (1 + t_k^2)), t_k = tan(pi/6 + pi k/14), on the real
cubic y^2 = x^2 (x - 1) have exactly the zero-sum triples of Z/14 as collinear triples, in Q(zeta_84).
"""
from itertools import combinations

N = 14

field = CyclotomicField(84)
zeta = field.gen()
i = zeta ** 21

points = []
for k in range(N):
    w = zeta ** (7 + 3 * k)
    t = (w - w ** -1) / (i * (w + w ** -1))
    assert t == t.conjugate()

    x = 1 + t * t
    y = t * x
    assert y * y == x * x * (x - 1)

    points.append((field(1), x, y))

assert len(set(points)) == N

triples = list(combinations(range(N), 3))
collinear = {T for T in triples if matrix(field, [points[p] for p in T]).det() == 0}
zero_sum = {T for T in triples if sum(T) % N == 0}

print(len(collinear), "collinear triples,", len(zero_sum), "zero-sum triples, equal:", collinear == zero_sum)
