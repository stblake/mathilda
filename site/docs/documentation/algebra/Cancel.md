# Cancel

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Cancel[expr] cancels out common factors in the numerator and denominator of expr.
Option Extension -> alpha cancels factors over Q(alpha) (e.g. simplifies
(x^2 - 2)/(x - Sqrt[2]) to x + Sqrt[2] when Extension -> Sqrt[2]).
Default Extension -> None treats algebraic numbers as opaque.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Cancel[(x^2 - 1) / (x - 1)]
Out[1]= 1 + x

In[2]:= Cancel[(x - y)/(x^2 - y^2) + (x^3 - 27)/(x^2 - 9)]
Out[2]= (9 + 3 x + x^2)/(3 + x) + 1/(x + y)

In[3]:= Cancel[(y - 1)/(Sqrt[y] - 1)]
Out[3]= 1 + Sqrt[y]

In[4]:= Cancel[(y - 1)/(y^(1/3) - 1)]
Out[4]= 1 + y^(1/3) + y^(2/3)

In[5]:= Cancel[1/(y^(2/3) - 1/y^(1/3))]
Out[5]= y^(1/3)/(-1 + y)

In[6]:= Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]]
Out[6]= Sqrt[2] + x

In[7]:= Cancel[(x^3 - 2)/(x - 2^(1/3)), Extension -> 2^(1/3)]
Out[7]= 2^(2/3) + 2^(1/3) x + x^2

In[8]:= Cancel[Sqrt[2]/(Sqrt[2] + Sqrt[2] x^4)]
Out[8]= 1/(1 + x^4)
```

## Implementation notes

**Algorithm.** `Cancel` reduces a rational expression to lowest terms by cancelling the polynomial GCD of numerator and denominator. `builtin_cancel_compute` strips an optional `Extension -> α` (with `Automatic` via `extension_autodetect`, routing single-generator and tower cases through the `Q(α)` paths `cancel_with_extension` / `qa_cancel_with_tower`) and otherwise calls `cancel_recursive`. That routine recurses through `List`/`Plus`/relational/logical heads, then for a leaf fraction: splits the expression into `num`/`den` with `extract_num_den` (which understands `Rational`, `Complex`, `Power`/`Exp` with negative or split exponents, and `Times`), strips common symbolic atoms (`rat_strip_symbolic_common`), computes `g = PolynomialGCD[num, den]`, and divides both sides by `g` using fraction-free exact polynomial division (`exact_poly_div` over the collected variable set). A soundness gate (`has_embedded_rational_subterm`) leaves the input untouched when a `Power[Plus/Times, negative]` is present, since the multivariate Euclidean GCD could blow up; the denominator's leading sign is normalised to positive.

**Data structures.** `Expr` trees throughout; the cancellation calls into the polynomial subsystem via the `PolynomialGCD` builtin (dispatching to `qaupoly_gcd` for extensions) and `exact_poly_div`, with variable sets collected by `collect_variables` and sorted by `compare_expr_ptrs`. `Numerator`/`Denominator` (same file) expose `extract_num_den` directly.

**Complexity / limits.** Dominated by the multivariate polynomial GCD; the embedded-rational gate and leaf-count fallbacks exist to avoid pathological Euclidean-blowup hangs, so some algebraically-cancellable inputs are returned uncancelled (left for later `Simplify` passes).

- `Protected`, `Listable`.
- Threads over equations, inequalities, logic functions, and sums dynamically.
- Evaluates greatest common divisors via polynomial GCD derivations avoiding extraneous expansions.
- Handles a single symbolic base appearing with rational fractional exponents (e.g. `Sqrt[y]`, `y^(1/3)`) by treating it as an algebraic generator: substitutes `y -> g^m` where `m` is the LCM of denominators, runs the polynomial cancellation in `g`, then substitutes back.
- The algebraic-generator pass runs `Together` on the substituted form (not just GCD-cancellation), so inputs whose `g`-substituted denominator is a Plus of terms with different `g`-denominators (e.g. `1/(g^2 - 1/g)` from `1/(y^(2/3) - 1/y^(1/3))`) are handled correctly.
- Extracts algebraic-constant atoms (`Sqrt[2]`, `Sqrt[3]`, `CubeRoot[5]`, `2^(2/3)`, ...) that appear in every summand of *both* numerator *and* denominator and divides them out before the polynomial GCD step. Closes a long-standing gap where `PolynomialGCD[Sqrt[2], Sqrt[2] + Sqrt[2] x^4]` returns 1 (the integer-content recursion treats `Sqrt[2]` as having content 1), so cancellations whose only shared factor was an algebraic constant survived as-is. The pass is intentionally narrow — only `Power[integer, rational/non-integer]` factors are eligible — to avoid disturbing the rational-function intermediates that the integration dispatcher pattern-matches against.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan) cancels common factors over `Q(alpha)` instead of `Q`. Implementation: lifts numerator and denominator into `Q(alpha)[x]` via the QAUPoly machinery, runs `qaupoly_gcd`, divides both sides by `g`, and re-renders. Works for single-fraction inputs; `Plus` inputs (sums of fractions) currently fall back to the no-extension path because `PolynomialQuotient` does not yet accept `Extension` (Phase 0.5 follow-up).

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- von zur Gathen & Gerhard, "Modern Computer Algebra", on polynomial GCD computation.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on rational function simplification.
- Source: [`src/rat.c`](https://github.com/stblake/mathilda/blob/main/src/rat.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Cancel[(x^2 - 1)/(x - 1)]
Out[1]= 1 + x
```

```mathematica
In[1]:= Cancel[(x^2 + 2 x + 1)/(x + 1)]
Out[1]= 1 + x
```

```mathematica
In[1]:= Cancel[(x^2 - 1)/(x^2 - 2 x + 1)]
Out[1]= (1 + x)/(-1 + x)
```

```mathematica
In[1]:= Cancel[(x^4 - 1)/(x^2 - 1)]
Out[1]= 1 + x^2
```

```mathematica
In[1]:= Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]]
Out[1]= Sqrt[2] + x
```

### Notes

Cancel removes the common polynomial factors between numerator and denominator by
dividing out their GCD, without performing a partial-fraction split. It reduces
`(x^2-1)/(x-1)` to `1 + x` and recognises perfect-square numerators such as
`(x^2+2x+1)/(x+1)`. When a common factor remains after cancellation, the result
stays a reduced quotient: `(x^2-1)/(x^2-2x+1)` becomes `(1+x)/(-1+x)` because both
share the factor `(x-1)` but the leftover `(x+1)/(x-1)` is already in lowest
terms. Cancel does not expand the surviving factors back out.

With the `Extension -> alpha` option, cancellation is performed over the
algebraic field `Q(alpha)`: `(x^2 - 2)/(x - Sqrt[2])` factors as
`(x - Sqrt[2])(x + Sqrt[2])` over `Q(Sqrt[2])`, so the denominator divides out
and the quotient collapses to `Sqrt[2] + x`. The default `Extension -> None`
treats algebraic numbers as opaque and would leave that quotient intact.
