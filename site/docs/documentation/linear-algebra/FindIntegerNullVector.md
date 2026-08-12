# FindIntegerNullVector

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FindIntegerNullVector[{x1, ..., xn}]`**

finds integers {a1, ..., an}, not all zero, with a1 x1 + ... + an xn == 0 (PSLQ / integer-relation detection).

**`FindIntegerNullVector[{x1, ..., xn}, d]`**

restricts the search to relations of norm \<= d.

<details>
<summary>Notes</summary>

The xi may be real or complex, exact or inexact; for complex xi the ai are Gaussian integers.  Exact relations are validated with PossibleZeroQ; for inexact xi the relation holds to the precision of the input.  When no relation is found the call is returned unevaluated. Options: WorkingPrecision    Automatic, or a digit count for the search. ZeroTest            Automatic, or a function applied to the residual.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= FindIntegerNullVector[{Log[2], Log[4]}]
Out[1]= {-2, 1}

In[2]:= FindIntegerNullVector[{Pi, ArcTan[1/5], ArcTan[1/239]}]
Out[2]= {1, -16, 4}

In[3]:= a = Sqrt[2] + 3^(1/3); FindIntegerNullVector[a^Range[0, 6]]
Out[3]= {1, -36, 12, -6, -6, 0, 1}

In[4]:= FindIntegerNullVector[{1, 2 I + Sqrt[3], (2 I + Sqrt[3])^2}]
Out[4]= {-7, -4*I, 1}
```

### Applications (4)

```mathematica
In[1]:= FindIntegerNullVector[{N[Zeta[2], 40], N[Pi^2, 40]}]
Out[1]= {-6, 1}
```

```mathematica
In[1]:= FindIntegerNullVector[{N[GoldenRatio, 40]^2, N[GoldenRatio, 40], 1}]
Out[1]= {-1, 1, 1}
```

```mathematica
In[1]:= FindIntegerNullVector[{N[Log[2], 40], N[Log[3], 40], N[Log[6], 40]}]
Out[1]= {-1, -1, 1}
```

```mathematica
In[1]:= FindIntegerNullVector[{N[Cos[Pi/7], 40]^3, N[Cos[Pi/7], 40]^2, N[Cos[Pi/7], 40], 1}]
Out[1]= {8, -4, -4, 1}
```

## Options & behaviour

> Implementation lives in `src/linalg/latticereduce.c` alongside
> `LatticeReduce`, reusing its exact Gaussian-rational LLL kernel
> (`lll_reduce`, extended to report `minᵢ |b*ᵢ|²`).

## Algorithm

latticereduce.c

LatticeReduce[{v1, v2, ...}] -- an LLL-reduced basis for the lattice spanned by the row vectors v_i.

```text
  LatticeReduce[m]   m an n x d matrix (List of n equal-length Lists).
                     Returns an n x d matrix whose rows form a reduced
                     basis of the same lattice (same Z- / Z[i]-module).
```

Entries may be:

```text
  - integers (machine int64 or GMP bigint),
  - rationals (Rational[p, q]),
  - Gaussian integers / Gaussian rationals (Complex[a, b] with a, b
    integer or rational).
```

Algorithm: the classical Lenstra-Lenstra-Lovas (LLL) reduction with

```text
Lovas parameter delta = 3/4, run in EXACT arithmetic.  Every quantity
```

is an exact Gaussian rational stored as a pair of GMP `mpq_t` (`GRat`), so the routine is correct for both machine-size and arbitrary-precision (bignum) inputs -- floating point is never used, which matters because LatticeReduce is most often used to discover integer relations where a rounding error would yield a wrong relation.

The Gram-Schmidt orthogonalisation is generalised to the complex

```text
(Hermitian) inner product  <x, y> = sum_k x_k conj(y_k),  so the same
code path handles real and Gaussian lattices.  Size reduction rounds
```

the mu coefficients to the nearest Gaussian integer; the Gram-Schmidt data (mu, |b*|^2) is maintained incrementally -- computed once up front, updated in place on size reduction, and updated on a Lovas swap via the conjugate-aware Cohen swap formulas (no full recomputation).

Because every basis transformation is an integer (Z, or Z[i] in the Gaussian case) row operation or a row swap, the lattice -- and hence Abs[Det] and every linear relation in the right null space -- is preserved exactly.

Linearly independent rows are required (every documented use of LatticeReduce -- integer-relation finding, basis reduction -- supplies

```text
a full-rank generating set).  A dependent generating set is detected
```

during Gram-Schmidt and the call is left unevaluated with a diagnostic.

```text
Memory ownership: standard builtin contract.  This file does NOT free
`res` on success or failure -- the evaluator owns it (see MEMORY.md /
```

SPEC.md S4.1).

## Implementation notes

**Algorithm.** `builtin_findintegernullvector` finds integers (or Gaussian integers, for complex inputs) `a = {a_1,…,a_n}`, not all zero, with `Σ a_i x_i = 0`, by **integer-relation detection via LLL** rather than PSLQ. It builds the relation lattice whose `i`-th row is `r_i = (e_i | round(2^b · x_i))` — the standard basis vector augmented with the scaled, rounded coordinate (Gaussian rounding when `x_i` is complex) — LLL-reduces it exactly using the same machinery as `LatticeReduce`, and reads the candidate relation off the leading components of the shortest reduced row. A rigorous certificate is computed from the LLL Gram–Schmidt bound `λ_1(L)² ≥ M2` (with `M2 = min_i |b*_i|²`) combined with the worst-case rounding error `|a·round(2^b x)| ≤ (√n/2)‖a‖`, giving a lower bound `B = √(M2 / (1 + (√n/2)²))` on the norm of any relation. `B` drives the no-relation / not-found diagnostics (`norel`, `lgrelb`, `rnfb`, `rnfu`).

**Data structures.** Exact Gaussian-rational `GRat` (pair of GMP `mpq_t`) scalars throughout, matching `LatticeReduce`. For inexact inputs the working precision `b` is derived from the inputs' Real/MPFR precision (`finv_max_prec_bits` under `USE_MPFR`).

**Complexity / limits.** Polynomial in `n` and the scaled bit-size; the exactness of the LLL pass plus the analytic norm bound let it both return a verified relation and *prove* none exists below a given norm (reported via the diagnostics).

- `Protected`.
- The `xi` may be **real or complex**, **exact or inexact**. For complex
  `xi` the `ai` are **Gaussian integers**.
- Built on the exact LLL machinery of `LatticeReduce`: the numbers are
  numericalised to a working precision `b` and embedded as the rows
  `(e_i | round(2^b x_i))` of an integer-relation lattice; the relation is
  read off the shortest reduced row.
- **Validation.** For exact `xi` the residual `a·x` is checked with
  `PossibleZeroQ`; a confidence guard (`n·log2(‖a‖²) < 1.35·b`) rejects
  large-coefficient artefacts that hold to only a few digits, forcing a
  precision increase. For inexact `xi` the relation holds to the precision
  of the input.
- **Precision** (`WorkingPrecision` option, default `Automatic`): inexact
  input uses the precision of the input; exact input starts at
  `$MachinePrecision` and escalates up to
  `$MachinePrecision + $MaxExtraPrecision` digits. An explicit digit count
  fixes the precision (no escalation).
- **Certified bound.** The rigorous LLL Gram–Schmidt bound
  `λ₁(L)² ≥ minᵢ |b*ᵢ|²` gives a lower bound
  `B = √(M² / (1 + n/4))` on the norm of *any* integer relation (`M² =
  minᵢ |b*ᵢ|²`), conditional on the numeric evaluation being correct to
  precision. `B` drives the no-relation diagnostics. Because this is the
  LLL bound (not the tighter PSLQ bound Wolfram reports), it is
  conservative: near the minimal-relation norm Mathilda may return the
  relation with a `lgrelb` message where Wolfram proves nonexistence.
- Relations are returned up to sign and (for complex input) up to a
  Gaussian-unit multiple — `{a}`, `{-a}`, and `{I a}` are all valid null
  vectors.
- Diagnostics: `FindIntegerNullVector::norel` (proven no relation with
  norm `≤ d`), `::lgrelb` (a relation was found but it exceeds `d`, and
  nonexistence is proven only below a smaller bound — the larger relation
  is returned), `::rnfb` (inexact/bounded: none found `≤ d`, proven none
  below a smaller bound), `::rnfu` (exact/unbounded: no relation found
  within the precision cap), `::ztest1` (the residual could not be
  proven zero and is assumed zero). When no vector is returned the call is
  left unevaluated.

**Attributes:** `Protected`.

## See also

[LatticeReduce](../../linear-algebra/LatticeReduce/), [PossibleZeroQ](../../expression-information/PossibleZeroQ/), [$MachinePrecision](../../expression-information/$MachinePrecision/)

## References

- A. K. Lenstra, H. W. Lenstra, L. Lovász, "Factoring Polynomials with Rational Coefficients", Mathematische Annalen 261 (1982).
- Henri Cohen, *A Course in Computational Algebraic Number Theory* (Springer, 1993).
- Source: [`src/linalg/latticereduce.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/latticereduce.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_findintegernullvector.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findintegernullvector.c)

## Notes & additional examples

### Notes

`FindIntegerNullVector` is integer-relation detection (PSLQ-style): given numerical values it recovers an exact integer combination summing to zero. Feeding `{ζ(2), π²}` recovers `−6 ζ(2) + π² = 0`, i.e. `ζ(2) = π²/6`. The powers of the golden ratio return `{-1, 1, 1}`, the minimal polynomial `−φ² + φ + 1 = 0`. The logarithms recover `log 6 = log 2 + log 3`. Most strikingly, the powers of `cos(π/7)` recover its minimal polynomial `8 x³ − 4 x² − 4 x + 1 = 0`, reconstructing exact algebraic structure from 40-digit floating-point samples. Working precision must comfortably exceed the size of the relation sought.
