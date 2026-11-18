import sys
import time
from itertools import combinations
from math import comb

from sage.all import GF, ZZ, matrix, prime_range


def zero_sum_lines(n):
    return [T for T in combinations(range(n), 3) if sum(T) % n == 0]


def boundary(n, lines):
    triples = {T: i for i, T in enumerate(combinations(range(n), 3))}
    quads = set()
    for line in lines:
        for d in range(n):
            if d not in line:
                quads.add(tuple(sorted(line + (d,))))
    entries = {}
    for r, Q in enumerate(sorted(quads)):
        for i in range(4):
            entries[(r, triples[Q[:i] + Q[i + 1:]])] = (-1) ** i
    return matrix(ZZ, len(quads), len(triples), entries, sparse=True)


def analyse(n, snf_limit):
    start = time.time()
    lines = zero_sum_lines(n)
    t = len(lines)
    M = boundary(n, lines)
    rank_q = M.change_ring(GF(1000003)).rank()
    free_defect = comb(n - 1, 3) - rank_q
    b3 = comb(n - 1, 2) - t
    phi3 = n * t - comb(n, 3) + b3 + free_defect
    torsion_primes = [p for p in prime_range(2, 2 * n + 2) if M.change_ring(GF(p)).rank() < rank_q]
    divisors = None
    if torsion_primes and M.nrows() <= snf_limit:
        divisors = [d for d in M.dense_matrix().elementary_divisors() if d > 1]
    print(f"n = {n}: t = {t}, phi_3 = {phi3}, decomposable bound 2t = {2 * t}, "
          f"Q-decomposable: {phi3 == 2 * t}, torsion primes {torsion_primes}, "
          f"torsion {divisors if divisors is not None else '-'} ({time.time() - start:.0f}s)", flush=True)


if __name__ == "__main__":
    low, high = int(sys.argv[1]), int(sys.argv[2])
    snf_limit = int(sys.argv[3]) if len(sys.argv) > 3 else 1500
    for n in range(low, high + 1):
        analyse(n, snf_limit)
