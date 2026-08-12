# PartitionsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PartitionsQ[n]`**

gives the number q(n) of partitions of the integer n into distinct

**`IntegerPartitions[n, All, Range[n]] with distinct parts.`**

<details>
<summary>Notes</summary>

parts (equivalently, into odd parts). n must be an integer; q(n) = 0 for n \< 0. Threads over lists. For the partitions themselves use

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Table[PartitionsQ[k], {k, 0, 20}]
Out[1]= {1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 27, 32, 38, 46, 54, 64}

In[2]:= PartitionsQ[100]
Out[2]= 444793

In[3]:= PartitionsQ[{2, 4, 6}]
Out[3]= {1, 2, 4}
```

## Algorithm

partitions.c — IntegerPartitions

A faithful, efficient recreation of the Wolfram-Language IntegerPartitions. The whole surface collapses onto a single count-vector enumerator over an ordered set of allowed parts, run with exact GMP rational arithmetic so that integers, big integers, rationals and negative values are all handled by the same code path.

Forms:

```text
  IntegerPartitions[n]                     all partitions of n
  IntegerPartitions[n, k]                  into at most k parts
  IntegerPartitions[n, {k}]                into exactly k parts
  IntegerPartitions[n, {kmin, kmax}]       between kmin and kmax parts
  IntegerPartitions[n, {kmin, kmax, dk}]   kmin, kmin+dk, ... parts
  IntegerPartitions[n, kspec, sspec]       parts drawn only from sspec
  IntegerPartitions[n, kspec, sspec, m]    first m (m>0) or last |m| (m<0)
```

n and the s_i may be rational and/or negative. Results are in reverse lexicographic order; within a partition the parts appear in the order of the reversed sspec (descending for the default Range[n]).

Ownership: this builtin only *reads* `res`. On every NULL return (bad arguments, ::undef, symbolic input) the evaluator keeps `res` unevaluated.

## Implementation notes

- `Protected`, `Listable` — `PartitionsQ[{2, 4, 6}]` → `{1, 2, 4}`.
- Two engines, dispatched by the size of `n` (threshold `n = 1000`):
  - **Small `n`** — an exact GMP recurrence derived from the Euler identity
    `∏(1−x^k)·∏(1+x^k) = ∏(1−x^{2k})`:
    `q(m) = r(m) + Σ_k (−1)^(k−1) [q(m − g_k) + q(m − g_k')]` over the same
    generalized pentagonal numbers `g_k = k(3k∓1)/2` as `PartitionsP`, with an
    inhomogeneous term `r(m) = (−1)^k` when `m/2` is the generalized pentagonal
    number `g_k` (else `0`). O(n^1.5) integer additions; O(n·√n)-bit table.
  - **Large `n`** — the **Hardy–Ramanujan–Rademacher / Hagis** exact convergent
    series evaluated in MPFR (Hagis, *Amer. J. Math.* 85 (1963) 213–222):
    `q(n) = (π/√(24n+1)) · Σ_{k odd} (1/k)·A_k(n)·I₁(π√(48n+2)/(12k))`, summing
    over **odd `k` only**, with the modified Bessel function `I₁` (computed by
    its everywhere-positive power series, as MPFR has no `I₁`) and the character
    sum `A_k(n) = Σ_{gcd(h,k)=1} cos(π(s(h,k) − s(2h,k)) − 2πnh/k)` built from
    exact Dedekind sums `s(h,k)`. Working precision `≈ π√(n/3)/ln 2` bits plus
    guard bits; memory stays tiny where the recurrence's table would be
    prohibitive.
- Unlike `PartitionsP`, the q-series has no clean closed-form Rademacher
  remainder bound, so the HRR engine's term count is convergence-driven (sum
  `≈√n` odd terms; accept only when the sum lies within `1/4` of an integer and
  the last term is below `1/8`, growing terms/precision otherwise). Correctness
  is pinned by an exhaustive HRR-equals-recurrence cross-check (every `n` up to
  1500 plus a strided sample and large spot values to `n = 100000`; see
  `tests/test_partitionsq.c`).
- `q(0) = 1`; `q(n) = 0` for `n < 0`. The result auto-promotes to a big integer
  (`q(1000)` has 22 digits; `q(100000)` has 245).
- Non-integer, symbolic, or astronomically large (big-integer) arguments are
  left unevaluated; `PartitionsQ` called with other than one argument emits
  `PartitionsQ::argx` and is left unevaluated.

**Attributes:** `Listable`, `Protected`.

## See also

[IntegerPartitions](../../number-theory/IntegerPartitions/), [PartitionsP](../../number-theory/PartitionsP/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_partitionsq.c`](https://github.com/stblake/mathilda/blob/main/tests/test_partitionsq.c)
