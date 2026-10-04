---
references:
  - "P. Hagis, *Partitions into odd summands*, Amer. J. Math. 85 (1963) 213-222."
  - "H. Rademacher, *Topics in Analytic Number Theory*, Springer, 1973."
source: src/partitions.c
---
**Algorithm.** `builtin_partitionsq` returns `q(n)`, the number of partitions of `n` into
distinct parts (equivalently, into odd parts; OEIS A000009), with the same size dispatch as
`PartitionsP` at a threshold of `n = 1000`. Below it, `partitionsq_recurrence` applies an
exact GMP recurrence derived from `prod(1 - x^k) Q(x) = prod(1 - x^(2k))`: the same
generalized-pentagonal homogeneous sum as `p(n)`, plus an inhomogeneous term `r(m) = (-1)^k`
whenever `m` is even and `m/2` is a generalized pentagonal number of index `k`. At or above
it, `partitionsq_hrr` evaluates the Hardy–Ramanujan–Rademacher / Hagis convergent series in
MPFR, `q(n) = (pi/sqrt(24n+1)) sum_(k odd) (1/k) A_k(n) I_1(pi sqrt(48n+2)/(12k))`, where
`A_k(n)` is a Dedekind-sum character sum and `I_1` is the modified Bessel function summed
from its everywhere-positive power series.

**Data structures.** The recurrence uses a table of `n+1` `mpz_t`; the HRR path works in
`mpfr_t` at `~pi sqrt(n/3)/ln2 + 96` bits, and the Dedekind sums `s(h, k)` are computed
exactly in `mpq_t` by `O(log k)` reciprocity. The result is an `Integer`/`BigInt`. There is
no ND/packed/`Compile` path, and `PartitionsQ` is `Listable`, threaded by the evaluator.

**Complexity / limits.** The recurrence is `O(n^1.5)` time and `O(n sqrt n)` bits. The Hagis
series has no clean closed-form remainder bound, so termination is convergence-driven: the
rounding is accepted only when the sum lies within `1/4` of an integer and the last term is
below `1/8` (so the unsummed tail cannot flip the rounding), growing the odd-`k` term count
and precision otherwise, with correctness pinned by an exhaustive HRR-equals-recurrence
cross-check in the unit tests. Big-integer arguments are rejected (unevaluated), `q(n) = 0`
for `n < 0`, and a wrong argument count emits `PartitionsQ::argx`.
