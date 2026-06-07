# Limit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Limit[f, x -> a]
    finds the limit of f as x approaches a.
Limit[f, {x1 -> a1, ..., xn -> an}]
    iterated limit, applied rightmost-first.
Limit[f, {x1, ..., xn} -> {a1, ..., an}]
    multivariate (joint) limit.
Limit[f, x -> a, Direction -> d]
    specifies the direction of approach:
      Reals or "TwoSided" -- default two-sided limit
      "FromAbove" or -1   -- approach from above (x -> a^+)
      "FromBelow" or +1   -- approach from below (x -> a^-)
      Complexes           -- limit over all complex directions

May return a finite value, Infinity, -Infinity, ComplexInfinity,
Indeterminate, Interval[{lo, hi}], or the original unevaluated
expression when the limit cannot be determined.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_limit` is a layered dispatcher; each layer either resolves the limit and short-circuits or hands the problem to the next. `builtin_limit` first rationalizes any inexact coefficients (the symbolic machinery is rational-coefficient only), mutes transient `Power::infy`/`Infinity::indet` warnings, then `builtin_limit_impl` normalizes the three calling forms — `Limit[f, x->a]`, `Limit[f, {x->a,…}]` (iterated), and `Limit[f, {x..}->{a..}]` (multivariate) — plus the `Direction`/`Assumptions` options and the `Direction -> I`/`Complexes` branch-cut post-pass.

The core `compute_limit` runs (in order): reciprocal-trig rewrite; at ±∞ a hyperbolic→exponential rewrite + `Expand`; an early bail on opaque/discontinuous heads (Floor, Ceiling, Sign, unknown `f[…x…]`); then the layer cascade — **Layer 1** structural fast paths (numeric-point substitution, then generic continuous substitution via `Together`); `Abs[g]` direction-aware rewrite; `ArcTan`/`ArcCot`-of-divergent; **Layer 3** rational-function `P(x)/Q(x)` shortcuts (leading-degree comparison); `Log` of a finite-limit inner; a Gruntz-lite dominant-summand `Log[sum]` reduction; `Log`+linear merge; term-wise `Plus` summation; **Layer 5.3** `f^g` logarithmic reduction (for exponents that depend on x, which `Series` can't expand); a bounded-envelope squeeze; **Layer 2** the `Series`-based workhorse (calls `Series[]` symbolically and reads the leading term); **Layer 5.1** L'Hospital's rule with leaf-count growth guardrails; **Layer 6** bounded-oscillation `Interval` returns; an atom-substitution recursion for `Power[b, e(x)]` subterms; and a last-resort one-sided-disagreement probe returning `Indeterminate`. Any layer returning NULL means "couldn't make progress here"; if all fail the `Limit[…]` is left unevaluated.

**Data structures.** A `LimitCtx { x, point, direction, depth }` threads the variable, the approach point, the (collapsed) direction, and a recursion counter (capped at `LIMIT_MAX_DEPTH`) through every layer. Subexpressions are manipulated as `Expr*` trees and most analysis is delegated to the evaluator (`Series`, `D`, `Together`, `Expand` are invoked symbolically through the symbol table, not via direct C calls).

**Complexity / limits.** `Series` subsumes most classical cases; L'Hospital is reserved for shapes `Series` can't expand and is guarded against complexity blow-up. Discontinuous-head and undefined-head inputs are deliberately refused rather than evaluated at a single side.

**Attributes:** `HoldAll`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 3.
- G. Gruntz, *On Computing Limits in a Symbolic Manipulation System*, PhD thesis, ETH Zürich, 1996.
- Source: [`src/calculus/limit.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/limit.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Limit[(Cos[x] - 1)/x^2, x -> 0]
Out[1]= -1/2
```

```mathematica
In[1]:= Limit[(1 + 1/x)^x, x -> Infinity]
Out[1]= E
```

```mathematica
In[1]:= Limit[(x^2 - 1)/(x - 1), x -> 1]
Out[1]= 2
```

```mathematica
In[1]:= Limit[Tan[x]/x, x -> 0]
Out[1]= 1
```

### Notes

`Limit[f, x -> a]` resolves the standard removable-singularity and indeterminate forms, including the classic `(1 + 1/x)^x -> E` and `0/0` cancellations such as `(x^2 - 1)/(x - 1)`. The `Direction` option selects one-sided (`"FromAbove"`/`"FromBelow"`) or complex approaches; the default is two-sided. Results may be a finite value, `Infinity`, `ComplexInfinity`, `Indeterminate`, an `Interval`, or the original expression unevaluated when the limit cannot be determined. Iterated and joint multivariate limits are supported through the list forms of the second argument.
