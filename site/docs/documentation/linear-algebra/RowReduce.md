# RowReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
RowReduce[m]
    gives the row-reduced form of the matrix m.
RowReduce[m, Method -> "<name>"]
    runs a specific elimination algorithm.  Accepted method names:
      "Automatic"                 — alias for "DivisionFreeRowReduction" (default)
      "DivisionFreeRowReduction"  — Bareiss-like fraction-free Gauss-Jordan
      "OneStepRowReduction"       — classical Gauss-Jordan with division per pivot
      "CofactorExpansion"         — identity-if-invertible via Laplace cofactor
                                     Det[m] (singular / rectangular m falls back
                                     to "DivisionFreeRowReduction")
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= RowReduce[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]
Out[1]= {{1, 0, -1}, {0, 1, 2}, {0, 0, 0}}

In[2]:= RowReduce[{{a, b, c}, {d, e, f}, {a+d, b+e, c+f}}]
Out[2]= {{1, 0, (c e)/(-b d + a e) - (b f)/(-b d + a e)}, {0, 1, -(c d)/(-b d + a e) + (a f)/(-b d + a e)}, {0, 0, 0}}

In[3]:= RowReduce[{{2, 1, 0}, {0, 3, 1}, {1, 0, 2}}, Method -> "CofactorExpansion"]
Out[3]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}

In[4]:= RowReduce[{{a, b}, {c, d}}, Method -> "OneStepRowReduction"]
Out[4]= {{1, 0}, {0, 1}}
```

## Implementation notes

**Algorithm.** `builtin_rowreduce` reduces `m` to reduced row-echelon form, dispatching on an optional `Method` option (`matsol_parse_method_option`):

- **`DivisionFreeRowReduction`** (the default / `Automatic`): `rowreduce_divfree`, a Bareiss-like fraction-free Gauss-Jordan. Per pivot column it picks the first structurally non-zero pivot (`is_zero_poly`), eliminates every other row with the Bareiss update `M[i,j] ← (pivot·M[i,j] − M[i,c]·M[r,j]) / P` where `P` is the previous pivot and the division is exact (`exact_div_wrapper`), keeping all intermediates polynomial/integer with no GCD blow-up. A final pass scales each pivot row to a leading 1, reducing each entry by its `PolynomialGCD` with the leading coefficient and `expand`-ing.
- **`OneStepRowReduction`**: `rowreduce_onestep`, textbook Gauss-Jordan — normalise the pivot row by dividing by the pivot (`matsol_div_entry`, which prefers exact division and canonicalises with `Together`), then subtract multiples from every other row. Result is already RREF.
- **`CofactorExpansion`**: `rowreduce_cofactor` returns identity-if-invertible for a non-singular square matrix, otherwise warns (`RowReduce::cofnsq`) and falls back to the fraction-free path.

**Data structures.** A flat row-major `Expr**` of size `m·n`; all arithmetic (`Times`, `Plus`, `Power`, `PolynomialGCD`, `Together`) is driven through the evaluator via `eval_and_free`, so integer, rational, complex, and free-symbolic matrices share one code path. Pivot detection is purely structural (`is_zero_poly`), so a `Real` cancellation of magnitude ~1e-18 is not treated as zero. Bad `Method` values emit `RowReduce::method` and leave the call unevaluated.

- `Protected`.
- Uses fraction-free division logic to perform exact algorithmic reduction across numerical, rational, and symbolics expressions natively avoiding division errors.
- Lives in `src/linalg/linsolve.c`; the helper primitives (Laplace cofactor determinant, exact polynomial division, tensor flatten / dimensions) are exposed from `src/linalg/util.c` and `src/linalg/det.c`.
- Accepts an optional `Method -> "<name>"` argument:
  - `Method -> Automatic` or `Method -> "Automatic"` (default) — alias for `"DivisionFreeRowReduction"`.
  - `Method -> "DivisionFreeRowReduction"` — Bareiss-like fraction-free Gauss-Jordan. Best for exact integer / rational / symbolic input — never produces a denominator larger than necessary.
  - `Method -> "OneStepRowReduction"` — classical Gauss-Jordan with one division per pivot per element. Each entry is canonicalised via `Together` so symbolic cancellations are still detected. Fast on numeric matrices.
  - `Method -> "CofactorExpansion"` — for a non-singular square matrix, returns the identity (verified via `Det[m] != 0` computed by Laplace cofactor expansion). On singular or rectangular input, falls back to `"DivisionFreeRowReduction"` and emits `RowReduce::cofnsq`.
- Unknown method names emit `RowReduce::method` and the call remains unevaluated.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- E. H. Bareiss, "Sylvester's Identity and Multistep Integer-Preserving Gaussian Elimination", Math. Comp. 22 (1968).
- Source: [`src/linalg/linsolve.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/linsolve.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
