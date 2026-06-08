# FindIntegerNullVector

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FindIntegerNullVector[{x1, ..., xn}]
    finds integers {a1, ..., an}, not all zero, with a1 x1 + ... + an xn == 0 (PSLQ / integer-relation detection).
FindIntegerNullVector[{x1, ..., xn}, d]
    restricts the search to relations of norm <= d.
The xi may be real or complex, exact or inexact; for complex xi the ai are Gaussian integers.  Exact relations are validated with PossibleZeroQ; for inexact xi the relation holds to the precision of the input.  When no relation is found the call is returned unevaluated.
Options:
    WorkingPrecision    Automatic, or a digit count for the search.
    ZeroTest            Automatic, or a function applied to the residual.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FindIntegerNullVector[{Log[2], Log[4]}]
Out[1]= {-2, 1}

In[2]:= FindIntegerNullVector[{Pi, ArcTan[1/5], ArcTan[1/239]}]
Out[2]= {-101842074, 1612412935, 397496960}

In[3]:= a = Sqrt[2] + 3^(1/3); FindIntegerNullVector[a^Range[0, 6]]
Out[3]= FindIntegerNullVector[{1, Sqrt[2] + 3^(1/3), (Sqrt[2] + 3^(1/3))^2, (Sqrt[2] + 3^(1/3))^3, (Sqrt[2] + 3^(1/3))^4, (Sqrt[2] + 3^(1/3))^5, (Sqrt[2] + 3^(1/3))^6}]

In[4]:= FindIntegerNullVector[{1, 2 I + Sqrt[3], (2 I + Sqrt[3])^2}]
Out[4]= {-7, -4*I, 1}

In[5]:= FindIntegerNullVector[{E, Pi}, 1000000]
Out[5]= FindIntegerNullVector[{E, Pi}, 1000000]
```

## Implementation notes

**Algorithm.** `builtin_findintegernullvector` finds integers (or Gaussian integers, for complex inputs) `a = {a_1,…,a_n}`, not all zero, with `Σ a_i x_i = 0`, by **integer-relation detection via LLL** rather than PSLQ. It builds the relation lattice whose `i`-th row is `r_i = (e_i | round(2^b · x_i))` — the standard basis vector augmented with the scaled, rounded coordinate (Gaussian rounding when `x_i` is complex) — LLL-reduces it exactly using the same machinery as `LatticeReduce`, and reads the candidate relation off the leading components of the shortest reduced row. A rigorous certificate is computed from the LLL Gram–Schmidt bound `λ_1(L)² ≥ M2` (with `M2 = min_i |b*_i|²`) combined with the worst-case rounding error `|a·round(2^b x)| ≤ (√n/2)‖a‖`, giving a lower bound `B = √(M2 / (1 + (√n/2)²))` on the norm of any relation. `B` drives the no-relation / not-found diagnostics (`norel`, `lgrelb`, `rnfb`, `rnfu`).

**Data structures.** Exact Gaussian-rational `GRat` (pair of GMP `mpq_t`) scalars throughout, matching `LatticeReduce`. For inexact inputs the working precision `b` is derived from the inputs' Real/MPFR precision (`finv_max_prec_bits` under `USE_MPFR`).

**Complexity / limits.** Polynomial in `n` and the scaled bit-size; the exactness of the LLL pass plus the analytic norm bound let it both return a verified relation and *prove* none exists below a given norm (reported via the diagnostics).

- `Protected`.
- The `xi` may be **real or complex**, **exact or inexact**. For complex

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- A. K. Lenstra, H. W. Lenstra, L. Lovász, "Factoring Polynomials with Rational Coefficients", Mathematische Annalen 261 (1982).
- Henri Cohen, *A Course in Computational Algebraic Number Theory* (Springer, 1993).
- Source: [`src/linalg/latticereduce.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/latticereduce.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
