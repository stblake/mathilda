# Limit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Limit[f, x -> a]`**

finds the limit of f as x approaches a.

**`Limit[f, {x1 -> a1, ..., xn -> an}]`**

iterated limit, applied rightmost-first.

**`Limit[f, {x1, ..., xn} -> {a1, ..., an}]`**

multivariate (joint) limit.

**`Limit[f, x -> a, Direction -> d]`**

specifies the direction of approach: Reals or "TwoSided" -- default two-sided limit "FromAbove" or -1   -- approach from above (x -\> a^+) "FromBelow" or +1   -- approach from below (x -\> a^-) Complexes           -- limit over all complex directions

**`Limit[f, x -> a, Method -> m]`**

selects the internal strategy: Automatic          -- (default) try all strategies in order "Substitution"     -- continuity, Abs kink, atom/one-sided probes "RationalFunction" -- degree comparison for P(x)/Q(x) "Series"           -- Taylor/Laurent/Puiseux leading term "LHospital"        -- L'Hospital's rule for 0/0 and Inf/Inf "Asymptotic"       -- dominant-term / log / exp reductions "Bounded"          -- squeeze and bounded-oscillation Interval "Oscillatory"      -- normal form c0 + Sum cj E^(I thetaj) at +-Inf "Gruntz"           -- Gruntz mrv algorithm for exp-log towers A named method leaves Limit unevaluated when it does not apply. Each method is also callable directly as Limit\`m\[f, x -\> a\].

**`Limit[f, x -> a, Assumptions -> assum]`**

uses the sign, magnitude, or domain of a symbolic parameter (also read from an ambient Assuming\[...\] / $Assumptions; the option wins) to decide otherwise-indeterminate parametric limits, e.g. Limit\[x^n, x -\> Infinity, Assumptions -\> n \> 0\] = Infinity and Limit\[a^x, x -\> Infinity, Assumptions -\> a \> 1\] = Infinity. With no informative assumption the result is the ordinary one.

<details>
<summary>Notes</summary>

May return a finite value, Infinity, -Infinity, ComplexInfinity, Indeterminate, Interval\[{lo, hi}\], or the original unevaluated expression when the limit cannot be determined.

</details>

## Examples (17)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0]
Out[1]= 1

In[2]:= Limit[(x^2 - 1)/(x - 1), x -> 1]
Out[2]= 2

In[3]:= Limit[(1 + 1/x)^x, x -> Infinity]
Out[3]= E

In[4]:= Limit[1/x, x -> 0]
Out[4]= ComplexInfinity

In[5]:= Limit[x^2 + y^2, {x, y} -> {1, 2}]
Out[5]= 5
```

### Options (3)

```mathematica
In[6]:= Limit[1/x, x -> 0, Direction -> "FromAbove"]
Out[6]= Infinity

In[7]:= Limit[Sin[x]/x, x -> 0, Method -> "Series"]
Out[7]= 1

In[8]:= Limit[(2 x^2 + 1)/(x^2 + x), x -> Infinity, Method -> "RationalFunction"]
Out[8]= 2
```

### Applications (9)

```mathematica
In[9]:= Limit[Sin[x]/x, x -> 0]
Out[9]= 1

In[10]:= Limit[(x^2 - 1)/(x - 1), x -> 1]
Out[10]= 2

In[11]:= Limit[(1 + a/x)^x, x -> Infinity]
Out[11]= E^a

In[12]:= Limit[(Sin[x] - x + x^3/6)/x^5, x -> 0]
Out[12]= 1/120

In[13]:= Limit[(x^x - x)/(1 - x + Log[x]), x -> 1]
Out[13]= -2

In[14]:= Limit[x - Sqrt[x^2 + x], x -> Infinity]
Out[14]= -1/2

In[15]:= Limit[x^2 + y^2, {x, y} -> {1, 2}]
Out[15]= 5

In[16]:= Limit[1/x, x -> 0, Direction -> "FromAbove"]
Out[16]= Infinity

In[17]:= Limit[1/x, x -> 0, Direction -> "FromBelow"]
Out[17]= -Infinity
```

## Algorithm

============================================================================ limit.c -- Symbolic limits for Mathilda. ============================================================================

Implements the Mathematica-style Limit built-in per the pipeline sketched in plans/limit_candidate_spec.md. The architecture is a layered dispatcher; each layer either resolves the limit and short-circuits or passes the problem down to the next layer:

```text
    Layer 0 -- Interface normalization (three calling forms, Direction).
    Layer 1 -- Cheap structural fast paths.
    Layer 2 -- Series-based evaluation (leverages Series[] natively).
    Layer 3 -- Rational-function dispatch (P(x)/Q(x) short-cuts).
    Layer 5 -- L'Hospital + logarithmic reduction heuristics.
    Layer 6 -- Bound analysis (Interval[] returns).
```

The series layer in Mathilda is powerful enough that it subsumes most classical DELIMITER cases. L'Hospital is reserved for those shapes where Series cannot compute a useful expansion (unknown heads, non- analytic inputs, etc.).

Memory conventions follow Mathilda standards: every helper that returns an Expr* returns a freshly-allocated tree owned by the caller. The top-level built-in returns a newly-allocated result on success (the evaluator frees the original `res` for us on a non-NULL return) or NULL to leave `res` unevaluated. In particular we never free `res` ourselves -- that would be a double-free against src/eval.c.

The module is intentionally layered with small single-purpose helpers so new test failures can be addressed by extending or swapping a single layer rather than re-plumbing the whole pipeline. =========================================================================

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Limit x (Log[x+1]-Log[x]) at Infinity | 1.24 s | 2.09 s | -- |
| Limit Sin[x]/x at 0 | 0.963 s | 1.37 s | -- |
| Limit (Exp[x]-1-x)/x^2 at 0 | 0.218 s | 1.63 s | -- |
| Series Exp[Sin[x]] to order 20 | -- | -- | -- |
| Series 1/(1-x-x^2) to order 60 | -- | -- | -- |
| Series Log[1+Sin[x]] to order 24 | -- | -- | -- |

## Implementation notes

**Algorithm.** `builtin_limit` is a layered dispatcher; each layer either resolves the limit and short-circuits or hands the problem to the next. `builtin_limit` first rationalizes any inexact coefficients (the symbolic machinery is rational-coefficient only), mutes transient `Power::infy`/`Infinity::indet` warnings, then `builtin_limit_impl` normalizes the three calling forms — `Limit[f, x->a]`, `Limit[f, {x->a,…}]` (iterated), and `Limit[f, {x..}->{a..}]` (multivariate) — plus the `Direction`/`Assumptions` options and the `Direction -> I`/`Complexes` branch-cut post-pass.

The core `compute_limit` runs (in order): reciprocal-trig rewrite; at ±∞ a hyperbolic→exponential rewrite + `Expand`; an early bail on opaque/discontinuous heads (Floor, Ceiling, Sign, unknown `f[…x…]`); then the layer cascade — **Layer 1** structural fast paths (numeric-point substitution, then generic continuous substitution via `Together`); `Abs[g]` direction-aware rewrite; `ArcTan`/`ArcCot`-of-divergent; **Layer 3** rational-function `P(x)/Q(x)` shortcuts (leading-degree comparison); `Log` of a finite-limit inner; a Gruntz-lite dominant-summand `Log[sum]` reduction; `Log`+linear merge; term-wise `Plus` summation; **Layer 5.3** `f^g` logarithmic reduction (for exponents that depend on x, which `Series` can't expand); a bounded-envelope squeeze; **Layer 2** the `Series`-based workhorse (calls `Series[]` symbolically and reads the leading term); **Layer 5.1** L'Hospital's rule with leaf-count growth guardrails; **Layer 6** bounded-oscillation `Interval` returns; an atom-substitution recursion for `Power[b, e(x)]` subterms; and a last-resort one-sided-disagreement probe returning `Indeterminate`. Any layer returning NULL means "couldn't make progress here"; if all fail the `Limit[…]` is left unevaluated.

**Method dispatch.** The `Method` option groups those layers into six named strategies — `"Substitution"`, `"RationalFunction"`, `"Series"`, `"LHospital"`, `"Asymptotic"` (the `ArcTan`/`Log`/dominant-term/`f^g` layers) and `"Bounded"` (envelope + oscillation) — plus `Automatic` (all of them). `parse_method` maps the option value to a `LIMIT_M_*` tag stored on `LimitCtx.method`. Inside `compute_limit` a `TRY(group, layer)` macro runs a layer only when it belongs to the selected group; the selection is enforced **only at the outermost call** (`method != Automatic && depth == 1`), so recursive sub-limits — one-sided probes, L'Hospital iterations, `Abs` splitting, the polar multivariate substitution — always inherit the full cascade. A named method that resolves nothing therefore leaves the whole `Limit` unevaluated, matching the "fail ⇒ unevaluated" contract; an unrecognised value emits `Limit::method` and is likewise left unevaluated.

**Data structures.** A `LimitCtx { x, point, direction, depth, method }` threads the variable, the approach point, the (collapsed) direction, a recursion counter (capped at `LIMIT_MAX_DEPTH`), and the selected `Method` tag through every layer. Subexpressions are manipulated as `Expr*` trees and most analysis is delegated to the evaluator (`Series`, `D`, `Together`, `Expand` are invoked symbolically through the symbol table, not via direct C calls).

**Complexity / limits.** `Series` subsumes most classical cases; L'Hospital is reserved for shapes `Series` can't expand and is guarded against complexity blow-up. Discontinuous-head and undefined-head inputs are deliberately refused rather than evaluated at a single side.

- `Protected`, `ReadProtected` (matches Mathematica; `Limit` does *not*
  hold its arguments, so `Limit[%, x -> Infinity]` sees the evaluated `f`).
- Options: `Direction -> Automatic`, `Assumptions -> Automatic`,
  `Method -> Automatic`.
- **Direction** selects the approach: `Reals`/`"TwoSided"` (default),
  `"FromAbove"` (or `-1`), `"FromBelow"` (or `+1`), or `Complexes`.
- **Method** restricts the internal strategy cascade to a single named
  group. `Automatic` runs the full cascade in order; a named method runs
  only that group and leaves `Limit` unevaluated if it does not apply
  (an unrecognised name emits `Limit::method` and is likewise left
  unevaluated). The restriction applies only to the outermost call —
  recursive sub-limits always use the full cascade. Methods:
  - `"Substitution"` — continuity / direct substitution, `Abs` kink
    resolution, atom-substitution and one-sided probes.
  - `"RationalFunction"` — leading-degree comparison for `P(x)/Q(x)`.
  - `"Series"` — Taylor / Laurent / Puiseux leading-term expansion.
  - `"LHospital"` — L'Hospital's rule with growth guardrails.
  - `"Asymptotic"` — dominant-term / `Log` / exponential reductions at
    infinity, including `f^g` via `Exp[g Log f]`, and the compose-at-infinity
    rule: for `f[g(x)]` whose inner argument diverges to `±Infinity`, apply the
    builtin's own value at Infinity (`Erf[Infinity] = 1`, `Tanh[Infinity] = 1`,
    `ArcTan[Infinity] = Pi/2`, `Gamma[Infinity] = Infinity`, …). Functions that
    do not self-evaluate there (oscillatory `Sin`, `Cos`) fall through and yield
    `Indeterminate`. Also splits a `Plus` at `±Infinity` whose term-wise sum
    fails because several summands each diverge: the finitely-converging terms
    are peeled off and the leftover group is re-limited together so its mutual
    divergences cancel — closing the real log-part of a rational antiderivative
    (`- b Log[1 - c x + x^2] + b Log[1 + c x + x^2] -> 0`), which is what lets
    `Integrate[1/(1 + x^4), {x, -Infinity, Infinity}]` close on the FTC path.
  - `"Bounded"` — squeeze envelope and bounded-oscillation `Interval`. Also
    covers a bounded base raised to a divergent positive power (`exp -> +Infinity`):
    with `B >= |base|` the pointwise magnitude bound, the limit is `0` when
    `Limit[B]` lies in `[0, 1)` (e.g. `(Sin[1/x]/2)^(1/x^2) -> 0` and the
    shrinking-bound `(x Sin[1/x]/2)^(1/x^2) -> 0` at `x -> 0`), and `Infinity`
    when the base is bounded below by a constant `> 1` and positive (e.g.
    `(2 + Sin[1/x]/2)^(1/x^2) -> Infinity`).
  - `"Oscillatory"` — the oscillatory normal form at `±Infinity`. `TrigToExp`
    plus `Expand` rewrites `f` as `c_0(x) + Sum_j c_j(x) E^(I theta_j(x))` with
    pairwise-distinct real phases carrying no constant term (a constant offset
    is folded into the amplitude, so `Cos[x]` and `Cos[x + 1]` share the phase
    `x`) and oscillation-free amplitudes. Distinct phases are asymptotically
    orthogonal, so the form decides the limit:
    - every `|c_j| -> 0` — the oscillation is squeezed away and the answer is
      `Limit[c_0]`: `Sin[x]/x -> 0`, `2 + Cos[x^2]/x -> 2`, and
      `Sin[x]^2 + Cos[x]^2 -> 1` where every oscillatory group cancels;
    - `Sum_j |c_j| / |c_0| -> r < 1` with `c_0 -> ±Infinity` — the oscillation
      cannot change the sign or the order of `f`, so the answer is `Limit[c_0]`:
      `x + Cos[x] -> Infinity`, `x^2 (2 + Cos[x]) -> Infinity` (`r = 1/2`);
    - one group strictly dominates all the others, its phase diverges and its
      amplitude has bounded argument — the intermediate value theorem produces
      two sequences with different limits, so **no limit exists** and the result
      is `Indeterminate`: `x Sin[x]`, `E^x Cos[x]`, `Sin[Log[x]]`;
    - every phase is a real polynomial of degree `>= 1` with numeric
      coefficients and `Limit[c_0]` is finite — then the Cesàro mean of `f` is
      `Limit[c_0]` and the Cesàro mean of `|f|^2` is `|Limit[c_0]|^2 +
      Sum_j (lim |c_j|)^2` (the cross terms die by van der Corput, every phase
      *difference* being a non-constant polynomial). A surviving oscillation
      therefore contradicts any finite limit, and `±Infinity` is excluded either
      because `f` is bounded or, for real `f`, by the window mean when
      `|c_j| = O(x^deg theta_j)`. Result `Indeterminate`:
      `Sin[x] + Cos[x]`, `Sin[x]^2`, `Cos[x] - Cos[x + 1]`, `Sin[x] Sin[x^2]`,
      and — the case with no dominant summand at all —
      `(Cos[x^2]/x^2 - Cos[(x+1)^2]/(x+1)^2) x^3`, which is asymptotically
      `2 x Sin[x^2 + x + 1/2] Sin[x + 1/2]`.

    A **finite** limit point reduces to this analysis through `x = a ± 1/t`
    with `t -> +Infinity` — an oscillation at a point is an oscillation at
    infinity in `t`, with the identical normal form — and a two-sided limit
    requires both sides to agree: `Sin[1/x]/x -> Indeterminate`,
    `Sin[1/x]^2 -> Indeterminate`, `x Sin[1/x] -> 0` at `x -> 0`.

    The layer abstains (leaving `Limit` unevaluated) whenever a hypothesis is
    not verifiable: a symbolic amplitude (`a Sin[x]`, since `a = 0` has the
    limit `0`), an amplitude that still carries an oscillation (`Tan`, `Sec` and
    `Csc` leave an exponential in a denominator), an envelope exactly equal to
    the smooth part (`x^2 (1 + Cos[x])`), or a phase that neither diverges nor
    is polynomial.
  - `"Gruntz"` — Gruntz's most-rapidly-varying (mrv) algorithm for exp-log
    functions (his 1996 ETH thesis). The whole function is expanded as a
    series in its most rapidly varying subexpression `w -> 0+`, which
    structurally avoids the intermediate expression swell that defeats
    bottom-up series methods on cancellation-heavy nested exponentials —
    e.g. `E^x (E^(1/x - E^-x) - E^(1/x)) -> -1`, `(3^x + 5^x)^(1/x) -> 5`,
    `x/Log[x^(Log[x]^(2/Log[x]))] -> Infinity`. It also runs as a last-resort
    fallback inside `Automatic` (after the series/L'Hospital layers), so these
    hard limits resolve without naming a method. Scope: real exp-log towers
    (`Exp`, `Log`, `Power`, `+`, `*`) plus tractable trig at a vanishing
    argument (`Sin`, `Cos`, …). A pre-processing pass also isolates the
    essential singularity of the semi-tractable special functions `Erf`,
    `Erfc`, `ExpIntegralEi`, `LogGamma`, `Gamma`, `PolyGamma[m, ·]`, `Zeta`,
    and the *modified* Bessel functions `BesselK[nu, ·]` / `BesselI[nu, ·]`
    (monotonic `Exp[∓z]` envelopes; the oscillatory `BesselJ`/`BesselY` are
    excluded) at infinity (thesis §5.2 transforms 5.6/5.7/5.8): each `F[g]` with
    `g -> ±oo` is replaced by its asymptotic expansion — obtained from
    `Series[F, {·, Infinity, n}]` — so the singular part becomes an explicit
    `Exp`/log head the mrv engine can handle, e.g.
    `ExpIntegralEi[x + E^-x] E^-x x -> 1` (thesis 5.4),
    `(Erf[x - E^-x] - Erf[x]) E^x E^(x^2) -> -2/Sqrt[Pi]`,
    `LogGamma[x]/(x Log[x]) -> 1`, `PolyGamma[x]/Log[x] -> 1`,
    `x PolyGamma[1, x] -> 1`, `x (PolyGamma[x] - Log[x]) -> -1/2` (digamma has a
    `Log[x]` growth head, `PolyGamma[m≥1]` decays like `x^-m`; the order `m`
    rides the 2-arg node fixed), `(Zeta[x] - 1) 2^x -> 1`,
    `Log[Zeta[x] - 1]/x -> -Log[2]` (`Zeta[x] = 1 + 2^-x + 3^-x + ...` collapses
    onto its exponential-scale Dirichlet head; a constant-exponent split
    normalises base-shifted terms like `2^-(x+1) -> (1/2) 2^-x` so same-class
    ratios `(Zeta[x]-1)/(Zeta[x+1]-1) -> 2` resolve), `Exp[x] Sqrt[x] BesselK[0, x]
    -> Sqrt[Pi/2]`, `BesselI[0, x] Exp[-x] Sqrt[x] -> 1/Sqrt[2 Pi]` (the Bessel
    order `nu` may be symbolic — the leading envelope is order-independent), and
    deep log-towers whose leading terms
    cancel, e.g. thesis 8.19
    `(Log[Log[x]+Log[Log[x]]]-Log[Log[x]])/Log[Log[x]+Log[Log[Log[x]]]] Log[x] -> 1`
    (each `Log` is expanded by factoring its `w`-pole out first, then the mrv
    `Series` runs in the positive log-scale `-Log[w]` so `Log[-Log[w]]` stays
    real). A separate pre-pass resolves `Max`/`Min` by *eventual dominance*:
    `Max[a, b]` eventually equals whichever argument is larger for large `x`,
    decided by the **leading-term** sign of `a - b` (so `Max[1/x, 2/x] = 2/x`
    though both `-> 0`), recursing so nested and factor-wrapped forms work —
    `Max[x, x^2] -> Infinity`, `x Max[1/x, 2/x] -> 2`, `Min[x, Log[x]] ->
    Infinity`, `Exp[x] Max[Exp[-x], Exp[-2x]] -> 1`. When the difference `a - b`
    involves a semi-tractable special function it is isolated to its asymptotic
    form before the sign is read, so `Max[PolyGamma[x], Log[x]] - Log[x] -> 0`
    and `Exp[x] Sqrt[x] Max[BesselK[0,x], BesselK[0,2x]] -> Sqrt[Pi/2]` resolve.
    A comparison that hinges on bounded oscillation (`Max[Sin[x], 2]`) has no
    leading-term sign and abstains.
    Still **not covered** (left unevaluated — never a wrong value): the
    thesis-8.31 `Gamma` Stirling difference (whose x^x-scale tower needs a
    deeper `Series` cancellation than the machinery reaches — a flagged `Series`
    limitation, not a wrong answer), `PolyGamma` with a symbolic order or at
    `-Infinity` (pole lattice), `Zeta` at `-Infinity` (trivial zeros), and —
    under an *exclusive*
    `Method -> "Gruntz"` — the oscillatory `BesselJ`/`BesselY` and bounded
    `Max`/`Min` cases below (the monotonic mrv engine cannot expand bounded
    oscillation). Under `Automatic` those resolve via a squeeze layer:
    `BesselJ[nu, x]`/`BesselY[nu, x] -> 0` (and as bounded factors, e.g.
    `BesselJ[0, x]/x -> 0`, while `x BesselJ[0, x]` correctly stays
    unevaluated — it has no limit), and `Max`/`Min` of a bounded oscillation
    against a dominating definite limit, e.g. `Max[Sin[x], 2] -> 2`,
    `Max[Sin[x], x] -> Infinity`, `Min[Cos[x], -x] -> -Infinity`.
- **Assumptions** carry the sign, magnitude, or domain of a symbolic
  parameter into the limit, so shapes that are indeterminate in the parameter
  resolve. The assumption is read from the `Assumptions -> ...` option, an
  enclosing `Assuming[...]`, or `$Assumptions` — the option wins over the
  ambient value — and is threaded through the whole cascade via the shared
  assumption engine (the same `AssumeCtx` `Simplify` and `PossibleZeroQ` use).
  With no informative assumption the result is byte-for-byte the
  assumption-free one. Reasoning is sound-only: a verdict is emitted just when
  the facts *entail* it; an unknown sign leaves `Limit` on its ordinary path.
  Decided shapes:
  - **exponent sign** — `Limit[x^n, x -> Infinity, Assumptions -> n > 0] =
    Infinity`, `n < 0 -> 0`; at `x -> 0` from above the two invert;
    `x -> -Infinity` resolves only for a known *even* exponent
    (`Element[n, Evens] && n > 0 -> Infinity`).
  - **real base magnitude** — `Limit[a^x, x -> Infinity, Assumptions -> a > 1]
    = Infinity`, `0 < a < 1 -> 0`, inverted at `x -> -Infinity`. A *real* base
    `> 1` gives a real `Infinity`, distinct from the
    `Assumptions -> Abs[a] > 1 -> ComplexInfinity` verdict (unknown phase).
  - **coefficient / log sign** — the leading coefficient of a `Plus` or
    rational, or `Sign[Log[b]]` for a symbolic base, is read from the facts:
    `Limit[c x, x -> Infinity, Assumptions -> c > 0] = Infinity` (`-Infinity`
    for `c < 0`), `Limit[x Log[a], x -> Infinity, Assumptions -> a > 1] =
    Infinity`, and the compose / one-sided-pole cases follow
    (`Limit[Tanh[c x], x -> Infinity, Assumptions -> c < 0] = -1`,
    `Limit[c/x, x -> 0, Direction -> "FromBelow", Assumptions -> c < 0] =
    Infinity`).
  - **growth ordering / divergent monomial** — a strict order fact ranks powers
    (`Limit[x^n/x^m, x -> Infinity, Assumptions -> n > m] = Infinity`), and a
    function of one divergent parametric monomial `x^a` (`a > 0`) is closed by
    the substitution `t = x^a` (`Limit[x^a/(x^a + 1), x -> Infinity,
    Assumptions -> a > 0] = 1`), which also feeds the bounded envelope
    (`Limit[Sin[x]/x^p, x -> Infinity, Assumptions -> p > 0] = 0`).
  - domain facts imply signs (`Element[n, PositiveIntegers]` gives `n > 0`), and
    the pre-existing `Assumptions -> Abs[B] R c` power dispatch (`|B| < 1 -> 0`,
    `|B| > 1 -> ComplexInfinity`, `|B| == 1 -> Indeterminate`) is unchanged.
- **Every method is also a head.** `Limit`Series[f, x -> a]` is exactly
  `Limit[f, x -> a, Method -> "Series"]`, so a strategy can be named without
  threading an option through — the natural way to ask "does *this* layer
  decide the limit?".  The heads are `Limit`Automatic`, `Limit`Substitution`,
  `Limit`RationalFunction`, `Limit`Asymptotic`, `Limit`Bounded`,
  `Limit`Series`, `Limit`LHospital`, `Limit`Gruntz` and `Limit`Oscillatory`,
  each with its own `Information` string and the attributes
  `{Protected, ReadProtected}`.  The two positional arguments and every other
  option (`Direction`, `Assumptions`) are forwarded untouched; a `Method`
  option is dropped, since the head already names the method.  An abstention
  echoes the head the user asked for
  (`Limit`RationalFunction[Sin[x]/x, x -> 0]` stays unevaluated) rather than
  falling back to the cascade.
- **Joint multivariate** limits at the origin or `+Infinity` are decided by a
  polar/spherical substitution: the integrand is `Simplify`-normalised in
  `r`/angle coordinates (cancelling common `r`-powers so buried `0/0` shapes
  like `ArcTan[y^2/(x^2 + x^3)]` do not fold to a spurious value), the radial
  `r`-limit is taken, and the resulting angular form is inspected — an
  angle-free constant is the limit, while two probe directions that disagree
  yield `Indeterminate`. Direction sampling of the original integrand is used
  only as a fallback when the polar analysis is inconclusive.
- May return a finite value, `Infinity`, `-Infinity`, `ComplexInfinity`,
  `Indeterminate`, an `Interval[{lo, hi}]`, or the original expression
  unevaluated when the limit cannot be determined.

**Attributes:** `Protected`, `ReadProtected`.

## References

**See also:** [Abs](../../arithmetic/Abs/), [Log](../../elementary-functions/Log/), [Sin](../../elementary-functions/Sin/), [Cos](../../elementary-functions/Cos/), [Plus](../../arithmetic/Plus/), [Interval](../../other-advanced/Interval/), [TrigToExp](../../elementary-functions/TrigToExp/), [Expand](../../algebra/Expand/)

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 3.
- G. Gruntz, *On Computing Limits in a Symbolic Manipulation System*, PhD thesis, ETH Zürich, 1996.
- Source: [`src/calculus/limit.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/limit.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
- Tests: [`tests/test_besselk.c`](https://github.com/stblake/mathilda/blob/main/tests/test_besselk.c)
- Tests: [`tests/test_gruntz.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz.c)
- Tests: [`tests/test_gruntz_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_gruntz_stress.c)
- Tests: [`tests/test_inexact_dispatch.c`](https://github.com/stblake/mathilda/blob/main/tests/test_inexact_dispatch.c)

## Notes & additional examples

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
