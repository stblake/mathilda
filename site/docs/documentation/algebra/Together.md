# Together

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Together[expr] combines fractions over a common denominator, then cancels.
Option Extension -> alpha runs the final cancellation over Q(alpha) so
algebraic-coefficient simplifications fire (e.g. 1/(x-Sqrt[2]) + 1/(x+Sqrt[2])
collapses to (2 x)/(x^2 - 2) under Extension -> Sqrt[2]).
Default Extension -> None combines and cancels over the rationals only.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Together[a/b + c/d]
Out[1]= (b c + a d)/(b d)

In[2]:= Together[x^2/(x^2 - 1) + x/(x^2 - 1)]
Out[2]= x/(-1 + x)

In[3]:= Together[1/x + 1/(x + 1) + 1/(x + 2) + 1/(x + 3)]
Out[3]= (6 + 22 x + 18 x^2 + 4 x^3)/(6 x + 11 x^2 + 6 x^3 + x^4)

In[4]:= Together[x^2/(x - y) - x y/(x - y)]
Out[4]= x

In[5]:= Together[y^(5/8)*(y^(19/8) - y^(73/24)/(y^(2/3) - 1/y^(1/3)))]
Out[5]= -y^3/(-1 + y)

In[6]:= Together[1/(x - Sqrt[2]) + 1/(x + Sqrt[2]), Extension -> Sqrt[2]]
Out[6]= (2 x)/(-2 + x^2)
```

## Implementation notes

**Algorithm.** `builtin_together` is a `FactorMemo`-keyed memoisation wrapper around `builtin_together_compute`, so repeated `Together[arg, opts]` calls issued during one `Simplify` (which pushes a `FactorMemo`) hit the cache; standalone calls see no overhead. The compute routine puts everything over a common denominator and cancels.

The core walker `together_recursive` recurses structurally. For a `Plus[t_1, ..., t_n]` it first `together`s each summand, splits each into numerator/denominator via `extract_num_den`, computes the iterated `PolynomialLCM` of the denominators, then for each term forms `lcm_den / den_i` using `cancel_exact_div_strict` (exact polynomial division in Q[vars]) — a *strict* combine: if any quotient is not exact (a sign that `PolynomialLCM` was only valid in some algebraic-extension ring, e.g. multi-radical sums), it bails cleanly and leaves the `Plus` uncombined rather than risk the `Power[Plus[...], -1]`-inside-numerator GCD blowup. When all quotients are exact, it sums the rescaled numerators, multiplies by `Power[lcm_den, -1]`, and runs `cancel_recursive` (PolynomialGCD-based cancellation) on the result. `List`/relational/logical heads thread through their args; other heads recurse and then cancel.

Crucially, `together_recursive` never expands `Power[Plus[...], n]` — `extract_num_den` only splits products, powers (pushing negative integer/rational exponents into the denominator), exponentials, complexes, and literal rationals; expansion is left to downstream poly-GCD with size gates. `builtin_together_compute` additionally handles `Extension -> α` (and `Extension -> Automatic`, which calls `extension_autodetect`): single-generator extensions route through `together_recursive_ext` over Q(α)[x]; multi-generator towers route through `qa_cancel_with_tower`, with nested-radical α prefiltered out and a "best generator" guided fallback (`pick_best_tower_generator`) gated on leaf count. A `poly_find_radical_gen` pass substitutes a radical generator to a fresh symbol, runs the no-extension path, and substitutes back. Inexact inputs are rationalised then re-numericalised (`internal_rationalize_then_numericalize`).

**Data structures.** `Expr*` trees; numerator/denominator pairs are plain `Expr*` out-parameters. Extension towers are `QATower`; exact division and LCM run on the multivariate-polynomial representation behind `exact_poly_div` / `PolynomialLCM` / `PolynomialGCD`. `Numerator` carries `ATTR_LISTABLE`.

**Complexity / limits.** Bounded by the polynomial GCD/LCM machinery on the combined fraction. The strict-quotient gate and the leaf-count gates are deliberate safeguards against multivariate-Euclid runaway on multi-radical / algebraic-extension inputs.

- `Protected`, `Listable`.
- Makes a sum of terms into a single rational function.
- Computes lowest common multiples (LCM) of denominators securely without unconditionally destroying pre-factored bases unnecessarily.
- Handles a single symbolic base appearing with rational fractional exponents (e.g. `y^(1/3)`, `y^(2/3)`, `y^(73/24)`) by treating it as an algebraic generator: substitutes `y -> g^m` where `m` is the LCM of denominators, runs the polynomial pipeline in `g`, then substitutes back.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan) combines into a single fraction with the standard combiner, then runs `Cancel[..., Extension -> alpha]` on the result so algebraic-coefficient cancellations fire. Effective on simple inputs like `1/(x - Sqrt[2]) + 1/(x + Sqrt[2])`, which collapses to `(2 x)/(x^2 - 2)`. Inputs whose summands themselves carry algebraic-coefficient denominators are deferred to Phase 0.5 (which will plumb `Extension` through `PolynomialLCM` / `PolynomialQuotient` / `together_recursive`).
- **Option `Extension -> Automatic` with polynomial radicands** (Phase E, 2026-05-25): when the input contains exactly one distinct radical `Sqrt[poly]` or `Power[poly, 1/q]` whose radicand is a polynomial in free symbols (e.g. `Sqrt[p+q]`, `Power[1+x^2, 1/3]`), `qa_cancel_with_poly_radical` substitutes the radical with a fresh symbol `S`, runs `Together`, reduces the numerator and denominator modulo `S^q - poly` via `PolynomialRemainder`, optionally rationalises via `PolynomialExtendedGCD` (when it shrinks the result), and substitutes back. Sample collapses: `Together[1/(x - Sqrt[p+q]) + 1/(x + Sqrt[p+q]), Extension -> Automatic]` returns `(2 x)/(-p - q + x^2)`; `Cancel[(x^2 - (p+q))/(x - Sqrt[p+q]), Extension -> Automatic]` returns `Sqrt[p+q] + x`. Multi-radical inputs (Cardano-style conjugate pairs) are rejected; the deeper limitation is documented in `docs/spec/changelog/2026-05-25.md`.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on rational function arithmetic and common denominators.
- von zur Gathen & Gerhard, "Modern Computer Algebra", on polynomial GCDs used in cancellation.
- Source: [`src/rat.c`](https://github.com/stblake/mathilda/blob/main/src/rat.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Together[1/x + 1/y]
Out[1]= (x + y)/(x y)
```

```mathematica
In[1]:= Together[1/(x-1) + 1/(x+1)]
Out[1]= (2 x)/(-1 + x^2)
```

```mathematica
In[1]:= Together[a/b + c/d]
Out[1]= (b c + a d)/(b d)
```

### Notes

Together combines a sum of fractions over a single common denominator, the inverse
operation of `Apart`. It computes the least common denominator and cancels common
polynomial factors via GCD, so `1/(x-1) + 1/(x+1)` collapses to `(2 x)/(x^2-1)`
with the cross terms cleared. Distinct symbolic denominators are simply multiplied,
giving `(b c + a d)/(b d)` for `a/b + c/d`. Unlike a raw expansion, Together keeps
the denominator factored where possible and does not multiply out `Power[Plus, n]`
factors unnecessarily.
