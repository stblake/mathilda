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
Limit[f, x -> a, Method -> m]
    selects the internal strategy:
      Automatic          -- (default) try all strategies in order
      "Substitution"     -- continuity, Abs kink, atom/one-sided probes
      "RationalFunction" -- degree comparison for P(x)/Q(x)
      "Series"           -- Taylor/Laurent/Puiseux leading term
      "LHospital"        -- L'Hospital's rule for 0/0 and Inf/Inf
      "Asymptotic"       -- dominant-term / log / exp reductions
      "Bounded"          -- squeeze and bounded-oscillation Interval
    A named method leaves Limit unevaluated when it does not apply.

May return a finite value, Infinity, -Infinity, ComplexInfinity,
Indeterminate, Interval[{lo, hi}], or the original unevaluated
expression when the limit cannot be determined.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0]
Out[1]= 1

In[2]:= Limit[(x^2 - 1)/(x - 1), x -> 1]
Out[2]= 2

In[3]:= Limit[(1 + 1/x)^x, x -> Infinity]
Out[3]= E

In[4]:= Limit[1/x, x -> 0]
Out[4]= ComplexInfinity

In[5]:= Limit[1/x, x -> 0, Direction -> "FromAbove"]
Out[5]= Infinity

In[6]:= Limit[x^2 + y^2, {x, y} -> {1, 2}]
Out[6]= 5

In[7]:= Limit[Sin[x]/x, x -> 0, Method -> "Series"]
Out[7]= 1

In[8]:= Limit[(2 x^2 + 1)/(x^2 + x), x -> Infinity, Method -> "RationalFunction"]
Out[8]= 2
```

## Implementation notes

**Algorithm.** `builtin_limit` is a layered dispatcher; each layer either resolves the limit and short-circuits or hands the problem to the next. `builtin_limit` first rationalizes any inexact coefficients (the symbolic machinery is rational-coefficient only), mutes transient `Power::infy`/`Infinity::indet` warnings, then `builtin_limit_impl` normalizes the three calling forms — `Limit[f, x->a]`, `Limit[f, {x->a,…}]` (iterated), and `Limit[f, {x..}->{a..}]` (multivariate) — plus the `Direction`/`Assumptions` options and the `Direction -> I`/`Complexes` branch-cut post-pass.

The core `compute_limit` runs (in order): reciprocal-trig rewrite; at ±∞ a hyperbolic→exponential rewrite + `Expand`; an early bail on opaque/discontinuous heads (Floor, Ceiling, Sign, unknown `f[…x…]`); then the layer cascade — **Layer 1** structural fast paths (numeric-point substitution, then generic continuous substitution via `Together`); `Abs[g]` direction-aware rewrite; `ArcTan`/`ArcCot`-of-divergent; **Layer 3** rational-function `P(x)/Q(x)` shortcuts (leading-degree comparison); `Log` of a finite-limit inner; a Gruntz-lite dominant-summand `Log[sum]` reduction; `Log`+linear merge; term-wise `Plus` summation; **Layer 5.3** `f^g` logarithmic reduction (for exponents that depend on x, which `Series` can't expand); a bounded-envelope squeeze; **Layer 2** the `Series`-based workhorse (calls `Series[]` symbolically and reads the leading term); **Layer 5.1** L'Hospital's rule with leaf-count growth guardrails; **Layer 6** bounded-oscillation `Interval` returns; an atom-substitution recursion for `Power[b, e(x)]` subterms; and a last-resort one-sided-disagreement probe returning `Indeterminate`. Any layer returning NULL means "couldn't make progress here"; if all fail the `Limit[…]` is left unevaluated.

**Method dispatch.** The `Method` option groups those layers into six named strategies — `"Substitution"`, `"RationalFunction"`, `"Series"`, `"LHospital"`, `"Asymptotic"` (the `ArcTan`/`Log`/dominant-term/`f^g` layers) and `"Bounded"` (envelope + oscillation) — plus `Automatic` (all of them). `parse_method` maps the option value to a `LIMIT_M_*` tag stored on `LimitCtx.method`. Inside `compute_limit` a `TRY(group, layer)` macro runs a layer only when it belongs to the selected group; the selection is enforced **only at the outermost call** (`method != Automatic && depth == 1`), so recursive sub-limits — one-sided probes, L'Hospital iterations, `Abs` splitting, the polar multivariate substitution — always inherit the full cascade. A named method that resolves nothing therefore leaves the whole `Limit` unevaluated, matching the "fail ⇒ unevaluated" contract; an unrecognised value emits `Limit::method` and is likewise left unevaluated.

**Data structures.** A `LimitCtx { x, point, direction, depth, method }` threads the variable, the approach point, the (collapsed) direction, a recursion counter (capped at `LIMIT_MAX_DEPTH`), and the selected `Method` tag through every layer. Subexpressions are manipulated as `Expr*` trees and most analysis is delegated to the evaluator (`Series`, `D`, `Together`, `Expand` are invoked symbolically through the symbol table, not via direct C calls).

**Complexity / limits.** `Series` subsumes most classical cases; L'Hospital is reserved for shapes `Series` can't expand and is guarded against complexity blow-up. Discontinuous-head and undefined-head inputs are deliberately refused rather than evaluated at a single side.

- `HoldAll`, `Protected`, `ReadProtected`.
- Options: `Direction -> Automatic`, `Assumptions -> Automatic`,

**Attributes:** `HoldAll`, `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 3.
- G. Gruntz, *On Computing Limits in a Symbolic Manipulation System*, PhD thesis, ETH Zürich, 1996.
- Source: [`src/calculus/limit.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/limit.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0]
Out[1]= 1
```

```mathematica
In[1]:= Limit[(x^2 - 1)/(x - 1), x -> 1]
Out[1]= 2
```

```mathematica
In[1]:= Limit[(1 + a/x)^x, x -> Infinity]
Out[1]= E^a
```

```mathematica
In[1]:= Limit[(Sin[x] - x + x^3/6)/x^5, x -> 0]
Out[1]= 1/120
```

```mathematica
In[1]:= Limit[(x^x - x)/(1 - x + Log[x]), x -> 1]
Out[1]= -2
```

```mathematica
In[1]:= Limit[x - Sqrt[x^2 + x], x -> Infinity]
Out[1]= -1/2
```

```mathematica
In[1]:= Limit[x^2 + y^2, {x, y} -> {1, 2}]
Out[1]= 5
```

```mathematica
In[1]:= Limit[1/x, x -> 0, Direction -> "FromAbove"]
Out[1]= Infinity

In[2]:= Limit[1/x, x -> 0, Direction -> "FromBelow"]
Out[2]= -Infinity
```

### Choosing a method

`Limit` evaluates by running a cascade of strategy layers in turn; each either
resolves the limit or hands the problem to the next. `Method -> m` restricts the
top-level call to a single strategy group. The default `Method -> Automatic`
runs the whole cascade in the order below. A named method computes *only* that
group; if it does not apply to the given expression, the `Limit` is left
unevaluated (an unrecognised method name is reported and also left unevaluated).

```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0, Method -> "Series"]
Out[1]= 1

In[2]:= Limit[(2 x^2 + 1)/(x^2 + x), x -> Infinity, Method -> "RationalFunction"]
Out[2]= 2

In[3]:= Limit[Sin[x]/x, x -> 0, Method -> "RationalFunction"]
Out[3]= Limit[Sin[x]/x, x -> 0, Method -> "RationalFunction"]
```

| `Method` | Strategy | Typical use |
|----------|----------|-------------|
| `Automatic` | run every strategy below, in order | default — best all-rounder |
| `"Substitution"` | continuity / direct substitution (via `Together`), `Abs` kink resolution, atom-substitution and one-sided probes | removable singularities, `Abs`, essential-singularity ratios |
| `"RationalFunction"` | leading-degree comparison for `P(x)/Q(x)` | rational functions at a point or at `Infinity` |
| `"Series"` | Taylor / Laurent / Puiseux expansion, reading the leading term | the workhorse — most `0/0` and `∞/∞` forms |
| `"LHospital"` | L'Hospital's rule with growth guardrails | `0/0`, `∞/∞` where `Series` cannot expand |
| `"Asymptotic"` | dominant-term / `Log` / exponential reductions at infinity, including `f^g` via `Exp[g Log f]` | limits at `Infinity`, `(1 + a/x)^x`, `Log`-of-sum |
| `"Bounded"` | squeeze / bounded-envelope to 0 and bounded-oscillation `Interval` returns | `Sin[x^2]/x`, oscillatory numerators |

The method restriction applies only to the outermost call: recursive
sub-limits — one-sided probes, L'Hospital iterations, `Abs` splitting — always
run the full cascade, so e.g. `Method -> "Series"` still resolves a two-sided
pole by falling back to its one-sided branches.

### Notes

`Limit[f, x -> a]` resolves the standard removable-singularity and indeterminate forms, including the classic `(1 + 1/x)^x -> E` and `0/0` cancellations such as `(x^2 - 1)/(x - 1)`. The `Direction` option selects one-sided (`"FromAbove"`/`"FromBelow"`) or complex approaches; the default is two-sided. The `Method` option (see above) selects a specific internal strategy, defaulting to `Automatic`. Results may be a finite value, `Infinity`, `ComplexInfinity`, `Indeterminate`, an `Interval`, or the original expression unevaluated when the limit cannot be determined. Iterated and joint multivariate limits are supported through the list forms of the second argument.
