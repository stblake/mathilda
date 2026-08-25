# PartitionsP

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PartitionsP[n]`**

gives the number p(n) of unrestricted partitions of the integer n.

<details>
<summary>Notes</summary>

n must be an integer; p(n) = 0 for n \< 0. Threads over lists. For the partitions themselves use IntegerPartitions\[n\].

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Table[PartitionsP[k], {k, 0, 12}]
Out[1]= {1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77}

In[2]:= PartitionsP[100]
Out[2]= 190569292

In[3]:= PartitionsP[4096]
Out[3]= 6927233917602120527467409170319882882996950147283323368445315320451

In[4]:= Table[Times @@ PartitionsP[Last /@ FactorInteger[n]], {n, 12}]
Out[4]= {1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 2}
```

### Worked examples (1)

```mathematica
In[5]:= PartitionsP[{2, 4, 6}]
Out[5]= {2, 5, 11}
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

- `Protected`, `Listable` — `PartitionsP[{2, 4, 6}]` → `{2, 5, 11}`.
- Two engines, dispatched by the size of `n` (threshold `n = 1000`):
  - **Small `n`** — Euler's pentagonal-number-theorem recurrence
    `p(m) = Σ_k (-1)^(k-1) [p(m - g_k) + p(m - g_k')]` over the generalized
    pentagonal numbers `g_k = k(3k∓1)/2`, using exact GMP integers only.
    O(n^1.5) integer additions; O(n·√n)-bit table.
  - **Large `n`** — the non-recursive **Hardy–Ramanujan–Rademacher** exact
    formula evaluated in MPFR (Johansson, arXiv:1205.5991): a convergent
    series `p(n) = Σ_{k=1}^N √(3/k)·(4/(24n−1))·A_k(n)·U(C/k)` with
    `U(x) = cosh x − sinh x / x`, `C = (π/6)√(24n−1)`. The exponential sums
    `A_k(n)` use Selberg's cosine-sum identity (no Dedekind sums). The number
    of terms `N` is chosen from Rademacher's rigorous remainder bound so the
    rounded sum is provably exact; working precision is `≈ C/ln 2` bits plus
    guard bits. Uses only O(√n)-bit precision, so memory stays tiny while the
    recurrence's table would be prohibitive.
- Both engines are cross-validated to agree exactly (see
  `tests/test_partitionsp.c`). The result auto-promotes to a big integer
  (`p(4096)` has 67 digits; `p(10000)` has 106).
- `p(0) = 1`; `p(n) = 0` for `n < 0`.
- Non-integer, symbolic, or astronomically large (big-integer) arguments are
  left unevaluated; `PartitionsP` called with other than one argument emits
  `PartitionsP::argx` and is left unevaluated.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [N](../../arithmetic/N/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_nsum.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nsum.c)
- Tests: [`tests/test_partitionsp.c`](https://github.com/stblake/mathilda/blob/main/tests/test_partitionsp.c)
