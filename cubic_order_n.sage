"""The points of order dividing n on the real cubic y^2 = x^2 (x - 1), as a line arrangement:
P_k = (1, 1 + t_k^2, t_k (1 + t_k^2)), t_k = tan(pi/6 + pi k/n).

  1. exact check in Q(zeta_{12n}): the points are real, lie on the cubic, and their collinear triples are exactly
     the zero-sum triples mod n;
  2. h_3 = gr_3(G) from the integral holonomy presentation, with the torsion certificate, as in h3_direct.py.

    sage cubic_order_n.sage N
"""
import sys
from itertools import combinations

n = int(sys.argv[1]) if len(sys.argv) > 1 else 14

order = 12 * n
field = CyclotomicField(order)
zeta = field.gen()
i = zeta ** (order // 4)

points = []
for k in range(n):
    w = zeta ** (n + 6 * k)
    t = (w - w ** -1) / (i * (w + w ** -1))
    assert t == t.conjugate()

    x = 1 + t * t
    y = t * x
    assert y * y == x * x * (x - 1)

    points.append((field(1), x, y))

assert len(set(points)) == n

lines = sorted({
    tuple(sorted((a, b, (-a - b) % n)))
    for a, b in combinations(range(n), 2)
    if (-a - b) % n not in (a, b)
})
collinear = {T for T in combinations(range(n), 3) if matrix(field, [points[p] for p in T]).det() == 0}
assert collinear == set(lines)

print(f"n = {n}: {n} real points on the cubic, {len(lines)} collinear triples = zero-sum triples (exact, Q(zeta_{order}))")

pairs = list(combinations(range(n), 2))
pair_index = {p: k for k, p in enumerate(pairs)}

hall_basis = []
for a, b, c in combinations(range(n), 3):
    hall_basis += [("b", a, b, c), ("c", a, b, c)]
for a, b in pairs:
    hall_basis += [("aab", a, b), ("bab", a, b)]
hall_index = {h: k for k, h in enumerate(hall_basis)}


def bracket(m, i, j):
    if i == j:
        return {}

    sign = 1
    if i > j:
        i, j, sign = j, i, -1

    if m == i:
        return {hall_index[("aab", i, j)]: sign}
    if m == j:
        return {hall_index[("bab", i, j)]: sign}

    a, b, c = sorted((m, i, j))
    if m == b:
        return {hall_index[("b", a, b, c)]: sign}
    if m == c:
        return {hall_index[("c", a, b, c)]: sign}
    return {hall_index[("b", a, b, c)]: sign, hall_index[("c", a, b, c)]: -sign}


ideal = matrix(ZZ, len(lines), len(pairs))
for row, (a, b, c) in enumerate(lines):
    ideal[row, pair_index[(b, c)]] += 1
    ideal[row, pair_index[(a, c)]] -= 1
    ideal[row, pair_index[(a, b)]] += 1

kernel = ideal.right_kernel_matrix()

entries = {}
row = 0
for r in range(kernel.nrows()):
    relation = {pairs[j]: int(kernel[r, j]) for j in range(len(pairs)) if kernel[r, j]}
    for m in range(n):
        result = {}
        for (a, b), coefficient in relation.items():
            for k, value in bracket(m, a, b).items():
                result[k] = result.get(k, 0) + coefficient * value

        for k, value in result.items():
            if value:
                entries[row, k] = value
        row += 1

matrix_h3 = matrix(ZZ, row, len(hall_basis), entries, sparse=False)
torsion = sorted(int(d) for d in matrix_h3.elementary_divisors() if d not in (0, 1))
rank = len(hall_basis) - matrix_h3.rank()

print(f"n = {n}: h_3 = gr_3(G) has free rank {rank} and torsion elementary divisors {torsion}")

cubes = list(combinations(range(n), 3))
cube_index = {cube: k for k, cube in enumerate(cubes)}

rows = []
for a, b, c in lines:
    for v in range(n):
        row = [0] * len(cubes)
        for (p, q), sign in (((b, c), 1), ((a, c), -1), ((a, b), 1)):
            if v in (p, q):
                continue
            order = -1 if ((v > p) + (v > q)) % 2 else 1
            row[cube_index[tuple(sorted((v, p, q)))]] += sign * order
        rows.append(row)

mu = matrix(ZZ, rows)
D, U, V = mu.smith_form()
V_inverse = V.inverse()

certified = []
for k in range(min(D.nrows(), D.ncols())):
    d = D[k, k]
    if d in (0, 1):
        continue

    z = V_inverse.row(k)
    u = U.row(k)
    f = V.column(k)
    assert u * mu == d * z
    for p in prime_divisors(d):
        assert all(x % p == 0 for x in mu * f)
        assert (z * f) % p == 1

    certified.append(int(d))

assert sorted(certified) == torsion
print(f"n = {n}: certificate in E^3/(E^1 I^2): torsion {sorted(certified)} verified")
