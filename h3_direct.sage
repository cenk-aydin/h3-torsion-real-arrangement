"""h_3 of the 14-point real arrangement, computed from the integral holonomy Lie algebra.

h = L(x_1, ..., x_n) / ideal(R), with R the annihilator of the degree-2 Orlik-Solomon ideal I^2 in Lambda^2.
Then h_3 = L_3 / [L_1, R], with L_3 in its Hall basis over Z. R is also built from the standard relations
[x_i, sum_{j in X} x_j] as a cross-check. H3_TEST=braid4 runs the braid arrangement A_3 instead.

The torsion is also certified in E^3 / (E^1 I^2), whose torsion equals that of h_3: for each torsion divisor d
an integral z with d z in E^1 I^2, and for each prime p | d a functional f with f(z) = 1 and f(E^1 I^2) = 0 mod p.

    sage h3_direct.sage
"""
import os
from itertools import combinations


def zero_sum_lines(n):
    return sorted({
        tuple(sorted((a, b, (-a - b) % n)))
        for a, b in combinations(range(n), 2)
        if (-a - b) % n not in (a, b)
    })


def braid_lines():
    edges = list(combinations(range(4), 2))
    index = {e: k for k, e in enumerate(edges)}
    lines = [
        tuple(sorted((index[(a, b)], index[(a, c)], index[(b, c)])))
        for a, b, c in combinations(range(4), 3)
    ]
    return len(edges), sorted(lines)


if os.environ.get("H3_TEST") == "braid4":
    n, lines = braid_lines()
else:
    n = 14
    lines = zero_sum_lines(n)

pairs = list(combinations(range(n), 2))
pair_index = {p: k for k, p in enumerate(pairs)}

hall_basis = []
for a, b, c in combinations(range(n), 3):
    hall_basis += [("b", a, b, c), ("c", a, b, c)]
for a, b in pairs:
    hall_basis += [("aab", a, b), ("bab", a, b)]
hall_index = {h: k for k, h in enumerate(hall_basis)}

assert len(hall_basis) == (n ** 3 - n) // 3


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


def h3(relations):
    entries = {}
    row = 0
    for relation in relations:
        for m in range(n):
            result = {}
            for (i, j), coefficient in relation.items():
                for k, value in bracket(m, i, j).items():
                    result[k] = result.get(k, 0) + coefficient * value

            for k, value in result.items():
                if value:
                    entries[row, k] = value
            row += 1

    matrix_h3 = matrix(ZZ, row, len(hall_basis), entries, sparse=True)
    divisors = matrix_h3.dense_matrix().elementary_divisors()
    torsion = sorted(int(d) for d in divisors if d not in (0, 1))
    return len(hall_basis) - matrix_h3.rank(), torsion


ideal = matrix(ZZ, len(lines), len(pairs))
for row, (a, b, c) in enumerate(lines):
    ideal[row, pair_index[(b, c)]] += 1
    ideal[row, pair_index[(a, c)]] -= 1
    ideal[row, pair_index[(a, b)]] += 1

kernel = ideal.right_kernel_matrix()
relations_annihilator = [
    {pairs[j]: int(kernel[i, j]) for j in range(len(pairs)) if kernel[i, j]}
    for i in range(kernel.nrows())
]

flats = [set(line) for line in lines]
flats += [set(p) for p in pairs if not any(set(p) <= set(line) for line in lines)]

relations_standard = []
for flat in flats:
    for i in flat:
        relation = {}
        for j in flat:
            if j == i:
                continue
            key, sign = ((i, j), 1) if i < j else ((j, i), -1)
            relation[key] = relation.get(key, 0) + sign
        relations_standard.append(relation)

standard = matrix(ZZ, len(relations_standard), len(pairs), {
    (k, pair_index[p]): value
    for k, relation in enumerate(relations_standard)
    for p, value in relation.items()
})
same_lattice = standard.row_module() == kernel.row_module()

def cubic_closure(n, lines):
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

    return matrix(ZZ, rows)


def certified_torsion(mu):
    D, U, V = mu.smith_form()
    V_inverse = V.inverse()

    torsion = []
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

        torsion.append(int(d))

    return sorted(torsion)


print(f"{len(lines)} lines; rank of R: {kernel.nrows()}; standard relations span the same lattice: {same_lattice}")
for name, relations in (("R = ann(I^2)", relations_annihilator), ("standard relations", relations_standard)):
    rank, torsion = h3(relations)
    print(f"h_3 with {name}: free rank {rank}, torsion elementary divisors {torsion}")

certified = certified_torsion(cubic_closure(n, lines))
assert certified == torsion
print(f"certificate in E^3/(E^1 I^2): torsion {certified} verified")
