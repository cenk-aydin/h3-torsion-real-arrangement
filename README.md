# Torsion in $\mathfrak h_3$ of a real line arrangement

[![verify](https://github.com/cenk-aydin/h3-torsion-real-arrangement/actions/workflows/verify.yml/badge.svg)](https://github.com/cenk-aydin/h3-torsion-real-arrangement/actions/workflows/verify.yml)

This repository contains a real line arrangement whose integral holonomy Lie algebra has torsion in degree 3.
Equivalently, the third lower central series quotient $\mathrm{gr}_3$ of its fundamental group is not torsion-free.

## The arrangement

On the real cubic $y^2 = x^2(x-1)$, take the 14 points of order dividing 14:

```math
P_k = \bigl(1,\ 1+t_k^2,\ t_k(1+t_k^2)\bigr), \qquad t_k = \tan\left(\frac{\pi}{6} + \frac{\pi k}{14}\right), \qquad k = 0, \dots, 13.
```

Three of them are collinear exactly when $k_1 + k_2 + k_3 \equiv 0 \pmod{14}$.
This gives 26 triple points and no point of higher multiplicity.
Read as planes in $\mathbb C^3$, the points define a real arrangement $\mathcal A$.

## Result

For the group $G$ of the complement of $\mathcal A$,

```math
\mathfrak h_3(\mathcal A) = \mathrm{gr}_3(G) \cong \mathbb Z^{65} \oplus \mathbb Z/2.
```

The points of order dividing 20 on the same cubic (57 triple points) give odd torsion as well:

```math
\mathfrak h_3 = \mathrm{gr}_3(G) \cong \mathbb Z^{153} \oplus \mathbb Z/2 \oplus \mathbb Z/6 \oplus \mathbb Z/6 \oplus \mathbb Z/30.
```

## The cubic family

The points of order dividing $n$ on the same cubic give a real arrangement of $n$ planes for every $n$, with $t$ triple points and no points of higher multiplicity.
The arrangement is $\mathbb Q$-decomposable exactly when $\varphi_3 = \operatorname{rank} \mathrm{gr}_3(G)$ equals $2t$.
None of the following is.

| $n$ | $t$ | $\varphi_3$ | $2t$ | torsion of $\mathrm{gr}_3(G)$ |
| ---: | ---: | ---: | ---: | --- |
| 14 | 26 | 65 | 52 | $\mathbb Z/2$ |
| 15 | 31 | 88 | 62 | $0$ |
| 16 | 35 | 89 | 70 | $(\mathbb Z/2)^2$ |
| 17 | 40 | 89 | 80 | $0$ |
| 18 | 46 | 154 | 92 | $0$ |
| 19 | 51 | 110 | 102 | $(\mathbb Z/2)^3 \oplus \mathbb Z/38$ |
| 20 | 57 | 153 | 114 | $\mathbb Z/2 \oplus (\mathbb Z/6)^2 \oplus \mathbb Z/30$ |
| 21 | 64 | 189 | 128 | $\mathbb Z/2$ |
| 22 | 70 | 182 | 140 | $\mathbb Z/3 \oplus \mathbb Z/(3 \cdot 11 \cdot 101)$ |
| 23 | 77 | 167 | 154 | $\mathbb Z/(23 \cdot 3631)$ |
| 24 | 85 | 298 | 170 | $\mathbb Z/3$ |
| 25 | 92 | 216 | 184 | $(\mathbb Z/5)^4 \oplus \mathbb Z/(3 \cdot 5 \cdot 7481)$ |
| 26 | 100 | 263 | 200 | $\mathbb Z/2 \oplus \mathbb Z/144 \oplus \mathbb Z/(2^4 \cdot 3^2 \cdot 13^2 \cdot 53 \cdot 1009)$ |
| 27 | 109 | 353 | 218 | $\mathbb Z/3 \oplus \mathbb Z/57$ |
| 28 | 117 | 324 | 234 | $(\mathbb Z/2)^8 \oplus \mathbb Z/6 \oplus (\mathbb Z/24)^2 \oplus \mathbb Z/168 \oplus \mathbb Z/(2^3 \cdot 3^2 \cdot 5 \cdot 7 \cdot 73 \cdot 109)$ |
| 29 | 126 | 277 | 252 | $(\mathbb Z/2)^2 \oplus \mathbb Z/(2 \cdot 13 \cdot 953) \oplus \mathbb Z/(2 \cdot 5 \cdot 13 \cdot 953 \cdot 8597 \cdot 4876873)$ |
| 30 | 136 | 499 | 272 | $(\mathbb Z/255)^2$ |

For $n \le 13$ there is no torsion.

## Running the code

The scripts need [SageMath](https://www.sagemath.org/) (tested with 10.6), and `is_decomposable.m2` needs [Macaulay2](https://macaulay2.com/) (tested with 1.24.11).

```bash
sage exact_real.sage
sage h3_direct.sage
sage cubic_order_n.sage 20
sage cubic_family.sage 14 30
M2 --script is_decomposable.m2 14 29
```

| Script | What it does |
| --- | --- |
| `exact_real.sage` | Checks in exact arithmetic that the 14 points are real points of the cubic and that their collinear triples are exactly the 26 zero-sum triples. |
| `h3_direct.sage` | Computes $\mathfrak h_3$ from the holonomy Lie algebra over $\mathbb Z$ and certifies the torsion. With `H3_TEST=braid4` it runs the braid arrangement $A_3$ instead (rank 10, no torsion). |
| `cubic_order_n.sage N` | Does all of the above for the points of order dividing $N$ (default 14). |
| `cubic_family.sage A B` | Computes $t$, $\varphi_3$, decomposability and the torsion of $\mathrm{gr}_3$ for $A \le n \le B$, from the cubic part of the integral quadratic closure, which has the same torsion. |
| `is_decomposable.m2 N Q` | Runs `isDecomposable` of the Macaulay2 package HyperplaneArrangements on the same matroid, realised over $\mathbb Z/Q$ with $Q \equiv 1 \bmod N$ prime. |

A GitHub Actions workflow runs all scripts with SageMath 10.6 on every push; the badge above links to the runs and their output.

## Method

The holonomy Lie algebra is $\mathfrak h = \mathbb L(x_1, \dots, x_n)/\langle R \rangle$, where $R \subset \Lambda^2$ is the annihilator of the degree-2 part $I^2$ of the Orlik–Solomon ideal.
In degree 3, $\mathfrak h_3 = \mathbb L_3/[\mathbb L_1, R]$, computed in the Hall basis over $\mathbb Z$ and reduced to Smith normal form.

The torsion is also certified independently in $E^3/(E^1 I^2)$, the cubic part of the quadratic closure of the integral Orlik–Solomon algebra, which has the same torsion as $\mathfrak h_3$.
For each torsion divisor $d$ the scripts construct an integral cubic $z$ with $dz \in E^1 I^2$ and, for each prime $p \mid d$, a functional that vanishes on $E^1 I^2$ modulo $p$ but not on $z$.
Both are checked by assertions.

## Use of AI

The search that led to these examples, the code and the computations were developed with the help of Claude (Anthropic), an AI assistant.
The results do not depend on trusting it: the scripts verify every claim above by exact computation and assertions.

## Background

Torsion in lower central series quotients of arrangement groups was known in higher degrees, for example for the MacLane arrangement in degrees 5 and 6, and for the real arrangements of Artal Bartolo, Guerville-Ballé and Viu-Sos in degree 4.
The examples here come from a search for real arrangements whose Orlik–Solomon algebra is quadratic over $\mathbb Q$ but not over $\mathbb F_2$.
They are not such examples: their quadratic closure has a free defect as well.
