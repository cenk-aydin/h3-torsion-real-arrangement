"""The points of order dividing n on y^2 = x^2 (x - 1): three are collinear iff their indices sum to 0 mod n.

    sage cubic_family.sage A B          README table rows for A <= n <= B
    sage cubic_family.sage N N --m2     one line in the format of is_decomposable.m2

M is the boundary map from the 4-sets containing a collinear triple to the triples. Its cokernel is the cubic part
of the integral quadratic closure up to a free summand, so
    tors gr_3(G) = tors coker M,
    d = C(n-1, 3) - rank M          (free rank of the defect of the quadratic closure),
    phi_3 = rank gr_3(G) = n t - C(n, 3) + b_3 + d,   b_3 = C(n-1, 2) - t.
The arrangement is Q-decomposable iff phi_3 = 2t, the sum of phi_3(F_2) = 2 over the t triple points.
"""
import sys
from itertools import combinations
from math import comb

from sage.all import GF, QQ, ZZ, factor, matrix, prime_range


def zero_sum_lines(n):
    return [T for T in combinations(range(n), 3) if sum(T) % n == 0]


def boundary(n, lines):
    triples = {T: i for i, T in enumerate(combinations(range(n), 3))}
    quads = sorted({tuple(sorted(line + (d,))) for line in lines for d in range(n) if d not in line})
    entries = {}
    for r, Q in enumerate(quads):
        for i in range(4):
            entries[(r, triples[Q[:i] + Q[i + 1:]])] = (-1) ** i
    return matrix(ZZ, len(quads), len(triples), entries, sparse=True)


def cyclic(d):
    if d < 1000:
        return f"\\mathbb Z/{d}"
    parts = [f"{p}^{e}" if e > 1 else f"{p}" for p, e in factor(d)]
    return "\\mathbb Z/(" + " \\cdot ".join(parts) + ")"


def torsion_tex(divisors):
    if not divisors:
        return "$0$"
    terms = []
    for d in sorted(set(divisors)):
        m = divisors.count(d)
        terms.append(f"({cyclic(d)})^{m}" if m > 1 else cyclic(d))
    return "$" + " \\oplus ".join(terms) + "$"


def family_member(n):
    lines = zero_sum_lines(n)
    t = len(lines)
    M = boundary(n, lines)
    rank = M.change_ring(QQ).rank()
    d = comb(n - 1, 3) - rank
    b3 = comb(n - 1, 2) - t
    phi3 = n * t - comb(n, 3) + b3 + d
    assert phi3 >= 2 * t

    elementary = M.dense_matrix().elementary_divisors()
    assert len([x for x in elementary if x != 0]) == rank
    divisors = [x for x in elementary if x > 1]
    for p in prime_range(2, 2 * n + 2):
        assert rank - M.change_ring(GF(p)).rank() == len([x for x in divisors if x % p == 0])
    return t, phi3, divisors


if __name__ == "__main__":
    low, high = int(sys.argv[1]), int(sys.argv[2])
    for n in range(low, high + 1):
        t, phi3, divisors = family_member(n)
        if "--m2" in sys.argv:
            print(f"n = {n}: phi_3 = {phi3}, local sum 2t = {2 * t}, isDecomposable: {str(phi3 == 2 * t).lower()}")
        else:
            print(f"| {n} | {t} | {phi3} | {2 * t} | {torsion_tex(divisors)} |", flush=True)
