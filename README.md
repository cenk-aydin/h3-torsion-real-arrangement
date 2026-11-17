# Torsion in $\mathfrak h_3$ of a real line arrangement

This repository contains a real line arrangement whose integral holonomy Lie algebra has torsion in degree 3.
Equivalently, the third lower central series quotient $\mathrm{gr}_3$ of its fundamental group is not torsion-free.
This answers Question 6.3 in A. Suciu, *On the topology and combinatorics of decomposable arrangements*,
[arXiv:2404.04784](https://arxiv.org/abs/2404.04784), negatively.

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

## Running the code

The scripts need [SageMath](https://www.sagemath.org/) (tested with 10.6).

```bash
sage exact_real.sage
sage h3_direct.sage
sage cubic_order_n.sage 20
```

| Script | What it does |
| --- | --- |
| `exact_real.sage` | Checks in exact arithmetic that the 14 points are real points of the cubic and that their collinear triples are exactly the 26 zero-sum triples. |
| `h3_direct.sage` | Computes $\mathfrak h_3$ from the holonomy Lie algebra over $\mathbb Z$ and certifies the torsion. With `H3_TEST=braid4` it runs the braid arrangement $A_3$ instead (rank 10, no torsion). |
| `cubic_order_n.sage N` | Does all of the above for the points of order dividing $N$ (default 14). |

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
