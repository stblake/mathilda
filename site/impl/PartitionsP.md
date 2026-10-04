---
references:
  - "F. Johansson, *Efficient implementation of the Hardy-Ramanujan-Rademacher formula*, LMS J. Comput. Math. 15 (2012) 341-359 (arXiv:1205.5991)."
  - "H. Rademacher, *On the partition function p(n)*, Proc. London Math. Soc. 43 (1937) 241-254."
source: src/partitions.c
---
**Algorithm.** `builtin_partitionsp` returns `p(n)`, the number of unrestricted partitions
of `n`, dispatched by size at a threshold of `n = 1000`. Below it, `partitionsp_recurrence`
uses Euler's pentagonal-number-theorem recurrence `p(m) = sum_k (-1)^(k+1) [p(m - g1_k) +
p(m - g2_k)]` over the generalized pentagonal numbers `g_{1,2} = k(3k∓1)/2`, in exact GMP
integers. At or above it, `partitionsp_hrr` evaluates the non-recursive
Hardy–Ramanujan–Rademacher exact formula in MPFR, summing `W_k(n) · U(C/k)` terms (with
`C = (pi/6) sqrt(24n-1)`) and growing the term count `N` until the Rademacher truncation
bound `M(n, N) < 1/4` guarantees that rounding the real sum to the nearest integer is exact,
retrying at doubled precision if the rounding margin is thin. Without MPFR the HRR engine is
unavailable and the recurrence is the fallback; `p(n) = 0` for `n < 0`.

**Data structures.** The recurrence allocates a table of `n+1` `mpz_t`; the HRR path works
in `mpfr_t` at `~C/ln2 + 64` bits, so its working memory is `O(sqrt n)` bits rather than the
recurrence's `O(n sqrt n)`. Johansson's Algorithm 1 advances the Rademacher residue with
integer adds only, so a cosine is evaluated solely at the `O(sqrt k)` solutions. The result
is an `Integer`/`BigInt`. There is no ND/packed/`Compile` path — the answer is a bignum —
and `PartitionsP` is `Listable`, threaded by the evaluator.

**Complexity / limits.** The recurrence is `O(n^1.5)` operations with `O(n sqrt n)` bits of
storage — the reason a large `n` switches to HRR, which is `O(n^(3/4))` cosine evaluations at
tiny memory. A big-integer argument is rejected (unevaluated) because no table of that size
can be built; symbolic or non-integer arguments stay unevaluated; a wrong argument count
emits `PartitionsP::argx`.
