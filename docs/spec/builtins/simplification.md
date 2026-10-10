# Simplification

Functions for transforming an expression into a simpler equivalent form, and
the assumption machinery they consult. The core entry point is `Simplify`;
`SimplifyCount` exposes its complexity measure, `Assuming` / `$Assumptions` /
`Element` supply assumptions, and `TransformationFunctions`, `ComplexityFunction`
and `$SimplifyDebug` tune or trace the search.

## Simplify
Performs a sequence of algebraic and other transformations on an expression and
returns the simplest form it finds.

- `Simplify[expr]` — searches the built-in transformation collection and returns
  the candidate of minimum complexity.
- `Simplify[expr, assum]` — simplifies using the assumptions `assum`.

`Simplify` runs a bounded heuristic search: it repeatedly applies a battery of
transforms to the expression and its subexpressions, scores every candidate with
a complexity measure, and keeps the lowest-scoring form. Transforms compose
across rounds, and a candidate that scores worse than its parent is pruned before
it can seed the next round.

**Features**:
- `Protected`. **Not** `Listable`: a `List` in the assumption position is a
  conjunction of facts (see below), not a threading axis.
- The built-in transformation collection tries `Together`, `Cancel`, `Expand`,
  `Factor`, `FactorSquareFree`, `Apart`, `TrigExpand`, `TrigFactor`, a
  `TrigToExp`/`ExpToTrig` roundtrip, and per-variable `Collect`, keeping the
  smallest result.
- The default complexity measure is `SimplifyCount` — total subexpression count
  plus the decimal-digit count of integer leaves — so `100 Log[2]` is preferred
  over its expanded `Log[2^100]` form.
- **`Root[...]` objects** are treated as the constant algebraic numbers they are:
  a qqbar pre-pass canonicalises constant-algebraic subexpressions
  (`Simplify[Root[1+#^4&,2]^3 + Root[1+#^4&,2]]` → `I Sqrt[2]`), and an expression
  whose coefficients are `Root` objects is simplified without the multivariate
  blow-up that treating each `Root` as a polynomial generator once caused. A
  rational-function identity with `Root` coefficients collapses to `0` when the
  `Root`s are radical-expressible (degree ≤ 4, or a binomial); e.g.
  `Simplify[D[Integrate[Sqrt[Tan[x]], x], x] - Sqrt[Tan[x]]]` → `0`.
- Threads manually over `List`, `Equal`, `Unequal`, `Less`, `LessEqual`,
  `Greater`, `GreaterEqual`, `And`, `Or`, and `Not`, carrying any options through
  into each sub-call.

```mathematica
In[1]:= Simplify[(x - 1)(x + 1)(x^2 + 1) + 1]
Out[1]= x^4

In[2]:= Simplify[3/(x + 3) + x/(x + 3)]
Out[2]= 1

In[3]:= Simplify[a x + b x + c]
Out[3]= c + (a + b) x

In[4]:= Simplify[Sin[x]^2 + Cos[x]^2]
Out[4]= 1

In[5]:= Simplify[2 Tan[x]/(1 + Tan[x]^2)]
Out[5]= Sin[2 x]

In[6]:= Simplify[(E^x - E^(-x))/Sinh[x]]
Out[6]= 2

In[7]:= Simplify[{Sin[x]^2 + Cos[x]^2, 3/(x + 3) + x/(x + 3)}]
Out[7]= {1, 1}
```

### Specialised transformations

Beyond the generic algebraic transforms, `Simplify` reaches a number of
canonical forms that the generic pipeline alone does not:

- **Cross-base radical fusion** — distinct positive-integer radicals sharing an
  exponent are combined inside a `Times`.
- **Roots of unity** — `(-1)^(p/q)` and `E^(I p Pi/q)` atoms are reduced modulo
  the relevant cyclotomic polynomial.
- **Exact trig at rational multiples of Pi** — `Sin`/`Cos`/`Tan`/`Cot`/`Sec`/`Csc`
  of `r Pi` (`r` rational) is an algebraic number (built from the root of unity
  `exp(I pi r)`), extending exact closed forms past the small-denominator tables.
  A variable-free product or combination folds to its exact value:
  `Cos[Pi/7] Cos[2 Pi/7] Cos[3 Pi/7] -> 1/8`,
  `Tan[Pi/7] Tan[2 Pi/7] Tan[3 Pi/7] -> Sqrt[7]`, and `RootReduce[Sin[Pi/7]]`
  is an exact `Root[...]`.
- **Radical denesting** — `Sqrt[A + Sqrt[B]]` and cube-root towers collapse via
  the half-sum identity when the result is cleaner.
- **Inverse trig / hyperbolic identities** — standard relations such as
  `Sin[ArcCos[x]] == Sqrt[1 - x^2]` reduce, and the complementary-angle pairs
  `ArcSin[x]+ArcCos[x]`, `ArcTan[x]+ArcCot[x]`, `ArcSec[x]+ArcCsc[x]` all fold to
  `Pi/2`.
- **Inverse-trig rational-angle combinations** — a `Z`-linear combination of
  `ArcTan`/`ArcSin`/`ArcCos` of constants plus a rational multiple of `Pi` that
  is identically zero is decided exactly: `exp(I e)` is built from the Euler
  closed forms and proven to equal `1` over the algebraic numbers (so
  `e == 0 mod 2 Pi`), with a numeric screen only selecting the `2 Pi` branch.
  Machin-type identities reduce — `4 ArcTan[1/5] - ArcTan[1/239] - Pi/4 -> 0`,
  `ArcSin[3/5] + ArcSin[5/13] - ArcSin[56/65] -> 0`. The hyperbolic constant
  combinations (`ArcSinh`/`ArcCosh`/`ArcTanh`) reduce likewise via `e^e == 1`.
- **Conditional identities under `Assumptions`** — identities true only on a
  region/domain reduce when the assumptions establish it, via general machinery
  (never ad-hoc), with a `Reduce`/CAD entailment bridge for the sign/domain
  conditions the lightweight provers cannot settle:
  - radical-product combine `Sqrt[a] Sqrt[b] -> Sqrt[a b]` when a base is
    provably non-negative (`Simplify[Sqrt[x-1] Sqrt[x+1] - Sqrt[x^2-1], x>1] -> 0`);
  - inverse-of-forward `ArcTanh[Tanh[x]] -> x` (`x in Reals`), the hyperbolic
    analogue of `ArcTan[Tan[x]] -> x` on `-Pi/2 < x < Pi/2`;
  - `Abs[z]^2` / `Re` / `Im` via `ComplexExpand` when the variables are real;
  - inverse-function **addition** identities on a region, e.g.
    `ArcTan[x]+ArcTan[y]-ArcTan[(x+y)/(1-x y)] -> 0` on `-1<x<1 && -1<y<1`, and the
    `ArcTanh`/`ArcSinh`/`ArcCosh` analogues — decided by a derivative-constancy
    engine (gradient `== 0`, branch-cut-free on the connected region by `Reduce`,
    value `0` at a sample point) or a pointwise forward-function recognizer.
    These are sound: nothing reduces without the region established (a wrong or
    too-large region, e.g. crossing `x y = 1`, is left unchanged).
- **Logarithm simplification** — `Log` of a positive rational is decomposed over
  its prime factors, and linear combinations of logs are fused
  (`Sum c_i Log[a_i] -> Log[Prod a_i^c_i]`).
- **Log-power symmetry under positivity** — `base^exp -> Exp[exp Log[base]]`
  when the base and every single-arg `Log` argument in the exponent are provably
  positive, exposing the symmetry `x^Log[y] = y^Log[x]` so it cancels
  (`Simplify[x^Log[y] - y^Log[x], x>0 && y>0] -> 0`). Kept only on a strict
  complexity win, so a standalone `x^Log[x]` is left unchanged; the two-sided
  positivity gate is required — a one-sided assumption (only `y>0`) correctly
  does **not** reduce, since the identity fails at `x<=0`.
- **Constant complex powers** — a variable-free `Power[c1, c2]` whose base is not
  a positive real folds to its principal value `Exp[c2 Log[c1]]`:
  `Simplify[I^I - E^(-Pi/2)] -> 0`, `Simplify[(-1)^I E^Pi] -> 1`. Taken only on a
  whole-expression complexity win, so `2^I`, standalone `I^I`, and surds such as
  `(-1)^(1/3)` are preserved.
- **Inverse-function logarithmic forms** — `ArcCos` joins `ArcSin`/`ArcSinh`/
  `ArcTanh` in the `TrigToExp`/`ExpToTrig` table, so
  `Log[x + I Sqrt[1-x^2]] = I ArcCos[x]` is recognised (principal-value general,
  unconditional). `ArcCosh[u] = Log[u + Sqrt[u^2-1]]` is recognised **under
  `u >= 1`** — the combined radical is branch-correct only there, so the gate is
  mandatory (it is left unreduced for `u <= -1`).
- **Range-gated inverse-of-direct trig** — `ArcSin[Sin[x]] -> x` on
  `-Pi/2 <= x <= Pi/2`, `ArcCos[Cos[x]] -> x` on `0 <= x <= Pi`,
  `ArcCot[Cot[x]] -> x` on `0 < x < Pi` (each inverse is a bijection from exactly
  that strip; mirrors the existing `ArcTan[Tan[x]]`). Off-range inputs are left
  unchanged.
- **Root-of-unity periodicity under integer assumptions** — for integer `n`,
  `E^(2 I Pi n + z) -> E^z`, `E^(I Pi n + z) -> (-1)^n E^z`, and
  `(-1)^(k + n) -> (-1)^k (-1)^n`, so e.g. `E^(x + 2 I Pi n) - E^x -> 0` and
  `Cos[n Pi] - E^(I n Pi) -> 0`.
- **Branch-gated power of an exponential** — `(E^w)^r -> E^(r w)` when `Im[w]` is
  provably in `(-Pi, Pi]` (so the principal `Log` recovers `w`), e.g.
  `Sqrt[E^(2 I x)] -> E^(I x)` under `-Pi/2 < x < Pi/2`. This is the sound,
  assumption-gated counterpart of `PowerExpand`'s unconditional collapse, so it
  does **not** reduce off the strip.
- **Pythagorean completion and reduction** for trig and hyperbolic squares.
- **Exact trig/exp zero-recognition** — a `Plus` that is a rational function of a
  single exponential kernel `t = E^(I x)` and is identically zero (canonically a
  Risch antiderivative diff-back `D[G] - f`) is proven `0` by exact rational
  point-evaluation on a Nullstellensatz grid — no numeric sampling, no slow
  `Together`. When that rigorous test declines — because of bare polynomial
  dependence on the kernel variable (e.g. the `x` in the `x E^x Sin[x]`
  diff-back) or mixed real+imaginary exponential kernels
  (`E^((1+I) x) = E^x E^(I x)`) — an exact `TrigToExp`-collapse fallback catches
  the identity: `Simplify[D[Integrate[x E^x Sin[x], x], x] - x E^x Sin[x]] -> 0`,
  and angle-addition identities such as
  `Sin[x] Cos[y] + Cos[x] Sin[y] - Sin[x + y] -> 0` collapse too. A circular trig
  head with an affine argument `k x + c Pi` (`c` rational) is handled by
  expanding the constant phase out of the kernel (the `c Pi` part contributes
  exact root-of-unity coefficients), so `Tan[Pi/2 - x] - Cot[x] -> 0`,
  `Tan[x] + Tan[x + Pi/3] + Tan[x + 2 Pi/3] - 3 Tan[3 x] -> 0`, and
  `Sin[x] Sin[Pi/3 - x] Sin[Pi/3 + x] - 1/4 Sin[3 x] -> 0` (and `PossibleZeroQ`
  agrees). The same applies to a hyperbolic head with an affine *imaginary*
  phase `k x + i c Pi` (`Cosh[i c Pi] = Cos[c Pi]`, `Sinh[i c Pi] = i Sin[c Pi]`
  fold to the same root-of-unity coefficients), so `Tanh[x + I Pi] - Tanh[x]
  -> 0`, `Tanh[I Pi/2 - x] + Coth[x] -> 0`, and `Tanh[x] + Tanh[x + I Pi/3] +
  Tanh[x + 2 I Pi/3] - 3 Tanh[3 x] -> 0`.
- **Trig / radical-trig rational normal form** — rational functions of trig and
  hyperbolic kernels are reduced to a canonical fraction modulo the Pythagorean
  ideal. A quadratic radical of a kernel (e.g. `Sqrt[Tan[x]]`, `Tan[x]^(3/2)`) is
  carried as an algebraic generator `l` with `l^2 = g`, so rational functions of
  `Sqrt[Tan[x]]` reduce too — `Simplify[Tan[x]/Sqrt[Tan[x]]] -> Sqrt[Tan[x]]`,
  and `D[Integrate[Sqrt[Tan[x]], x], x] // Simplify -> Sqrt[Tan[x]]`. Radicands
  that are rational with an *odd* generator in the denominator (`Cot = Cos/Sin`,
  `Csc = 1/Sin`) are handled too: the inverse odd-generator powers the relation
  injects are cleared before the denominator is rationalised, so
  `Simplify[Cot[x]/Sqrt[Cot[x]]] -> Sqrt[Cot[x]]`.
- **Multi-generator radical rational normal form** — a rational function of two
  or more distinct radical bases (e.g. `a^(1/3)` and `(a+b x)^(1/3)`) is reduced
  in the quotient ring `K[g_1, ..., g_n] / <g_k^{q_k} - base_k>`: each base is
  carried as an algebraic generator, the terms are combined over a common
  denominator, reduced modulo the generator relations, and the denominator is
  rationalised, before the radicals are substituted back. This recovers cross-base
  cancellations the single-generator `Together`/`Cancel` path cannot — e.g.
  `D[Integrate[1/(x^3 (a+b x)^(1/3)), x], x] // Simplify -> 1/(x^3 (a+b x)^(1/3))`.
  **Algebraic *constant* radicals** such as `Sqrt[2]` are collected as generators
  too (with the relation `s^q - c`, e.g. `s^2 - 2`), so a rational function of `x`,
  one `x`-dependent radical, and a constant radical in its coefficients — the shape
  of an antiderivative re-differentiated against its integrand — reduces too:
  `D[Integrate[(x^2+1)/(x^3 Sqrt[2x^4-2x^2+1]), x], x] - integrand // Simplify ->
  0` (previously this fell into the algebraic-field `Together` over
  `Q(x)[R]/(R^2-q)` and did not terminate). The pass engages only when at least
  one base carries a free symbol — a purely numeric radical identity is left to the
  `RootReduce`/qqbar pass — and the result is adopted only when its `SimplifyCount`
  strictly improves, so it never regresses a case.
- **Equation / inequality rebalancing** — a binary relation is normalised by
  dividing through the GCD of integer coefficients and partitioning terms across
  the relation; the rebalanced form is kept when its `SimplifyCount` is lower.
  Strict inequalities flip when divided by a negative.

```mathematica
In[1]:= Simplify[Sqrt[2] Sqrt[3]]
Out[1]= Sqrt[6]

In[2]:= Simplify[Sqrt[6] - Sqrt[2] Sqrt[3]]
Out[2]= 0

In[3]:= Simplify[1 - (-1)^(1/3) + (-1)^(2/3)]
Out[3]= 0

In[4]:= Simplify[1 - (-1)^(1/5) + (-1)^(2/5) - (-1)^(3/5) + (-1)^(4/5)]
Out[4]= 0

In[5]:= Simplify[Sqrt[3 + 2 Sqrt[2]] - (1 + Sqrt[2])]
Out[5]= 0

In[6]:= Simplify[Sin[ArcCos[x]] - Sqrt[1 - x^2]]
Out[6]= 0

In[7]:= Simplify[ArcSin[x] + ArcCos[x] - Pi/2]
Out[7]= 0

In[8]:= Simplify[Log[72] - 3 Log[2] - 2 Log[3]]
Out[8]= 0

In[9]:= Simplify[4 Sin[x]^2 Cos[x]^2 + 4 Sin[x] Cos[x] + 1]
Out[9]= 1/2 (3 - Cos[4 x] + 4 Sin[2 x])

In[10]:= Simplify[2 x - 4 y + 6 z - 10 == -8]
Out[10]= x + 3 z == 1 + 2 y

In[11]:= Simplify[-2 x < 4]
Out[11]= x > -2
```

### Assumptions

`Simplify[expr, assum]` simplifies under `assum`, which may be equations,
inequalities, domain specifications such as `Element[x, Integers]`, or logical
combinations of these. A list `{a1, a2, ...}` is treated as the conjunction
`And[a1, a2, ...]`. Under provable positivity / reality, `Simplify` applies
`Log`/`Power` identities — `Log[a b] -> Log[a] + Log[b]`, `(a b)^c -> a^c b^c`,
`(a^p)^q -> a^(p q)`, `Log[a^p] -> p Log[a]` and the like — whenever the
operand-domain conditions follow from the assumption set. Per-symbol sign facts
drive `Sqrt[x^2] -> x` / `-x` / `Abs[x]`, and integer facts drive the
`Sin[n Pi] -> 0`, `Cos[n Pi] -> (-1)^n` family.

A **deep sign oracle** extends `Abs[g] -> ±g` and `Sqrt[g^2] -> ±g` to an
*arbitrary* real `g` whose sign is fixed on the assumed region — not just bare
symbols. It layers the structural provers, a decomposition of the pole-bearing
trig family (`sign[Tan] = sign[Sin] sign[Cos]`, `sign[Sec] = sign[Cos]`, …), and
a sound Reduce/CAD base case (a sign is proved only when the opposite strict
inequality is unsatisfiable on the region). It declines — leaving `Abs`/`Sqrt`
intact — wherever the sign is not constant (e.g. `Abs[Cos[x]]` on `0 < x < Pi`,
which flips at `Pi/2`, or an unbounded region).

```mathematica
In[1]:= Simplify[Abs[Sin[x]], 0 < x < Pi/2]
Out[1]= Sin[x]

In[2]:= Simplify[Abs[Tan[x]] + Tan[x], Pi/2 < x < Pi]
Out[2]= 0
```

When no positional assumption and no `Assumptions` option are given, `Simplify`
reads the current value of `$Assumptions`.

```mathematica
In[1]:= Simplify[Sqrt[x^2], x > 0]
Out[1]= x

In[2]:= Simplify[Sqrt[x^2], Element[x, Reals]]
Out[2]= Abs[x]

In[3]:= Simplify[Log[a b], a > 0 && b > 0]
Out[3]= Log[a] + Log[b]

In[4]:= Simplify[(a^p)^q, a > 0 && Element[p, Reals]]
Out[4]= a^(p q)

In[5]:= Simplify[Sqrt[x^2 y^2], x > 0 && y < 0]
Out[5]= -x y

In[6]:= Simplify[Cos[k Pi], Element[k, Integers]]
Out[6]= (-1)^k

In[7]:= Simplify[Log[E^(x + y)], {Element[x, Reals], Element[y, Reals]}]
Out[7]= x + y
```

A predicate that appears literally among the assumed facts folds to `True`:

```mathematica
In[8]:= Simplify[x > 0, x > 0]
Out[8]= True
```

### Options

- **`Assumptions`** (default `$Assumptions`) — the facts assumed while
  simplifying. An explicit `Assumptions -> X` overrides `$Assumptions`; a
  positional assumption is conjoined with `$Assumptions`.
- **`ComplexityFunction`** (default the built-in `SimplifyCount` measure) — ranks
  candidate forms; `Simplify` returns the lowest-scoring one. A custom function
  `f` must return an integer or bigint for `f[candidate]`; otherwise the default
  is used. `ComplexityFunction -> Automatic` is a synonym for the default and
  takes the fast native scoring path. Compared with the default, `LeafCount`
  drops the integer-digit penalty.
- **`TransformationFunctions`** (default `Automatic`) — the functions applied to
  try to transform parts of `expr` (see [TransformationFunctions](#transformationfunctions)).
- **`TimeConstraint`** (default `Infinity`) — a wall-clock budget in seconds.
  A scalar `t` bounds the **whole call**: the budget is honoured on every path —
  the heuristic search, the specialised rational/polynomial/log-exp pipelines
  (including the SHAPE_RATIONAL input that dispatches straight to a pipeline),
  the bottom-up descent, and the seed phase — so the entire `Simplify[expr, …]`
  returns the best form found so far once `t` is exhausted. A list
  `{tLoc, tTot}` additionally caps each individual sub-expression's search at
  `tLoc` while the whole call is capped at `tTot` (matching `FullSimplify`); a
  one-element `{tLoc}` sets the per-sub-expression budget only. The check is
  synchronous between search steps, so it fails gracefully with no memory leak,
  but it does **not** interrupt a single long-running kernel call mid-flight (a
  pathological single `Together`/`Factor` is instead bounded by the poly engine's
  own degree/size guards). `Infinity` (the default) imposes no limit.

```mathematica
In[1]:= Simplify[1/(x - 1) + 1/(1 - x), TransformationFunctions -> {Together}]
Out[1]= 0

In[2]:= Simplify[Sin[x]^2 + Cos[x]^2, TransformationFunctions -> {Together}]
Out[2]= Cos[x]^2 + Sin[x]^2

In[3]:= Simplify[Sin[x]^2 + Cos[x]^2, TransformationFunctions -> {Automatic, Together}]
Out[3]= 1

In[4]:= Simplify[a + b, TransformationFunctions -> {(# /. a -> 0 &)}]
Out[4]= b
```

## FullSimplify
Tries a wider range of transformations than `Simplify`, drawing on a library of
known elementary- and special-function identities, and returns the simplest form
it finds.

- `FullSimplify[expr]` — simplify using the relevant function identities.
- `FullSimplify[expr, assum]` — simplify under the assumptions `assum`.

`FullSimplify` is a wrapper around `Simplify`: it scans `expr` for the function
heads it contains, gathers the transformation functions registered for those
heads, and feeds them to `Simplify` via `TransformationFunctions ->
{Automatic, ...}`. Because `Simplify` only ever keeps a candidate of strictly
lower complexity, **`FullSimplify` always yields at least as simple a form as
`Simplify`** but may take longer.

Identities are organised into per-function-family libraries that are loaded
**only when a matching head appears** in the input, so a request involving no
special functions costs no more than `Simplify`, and the collection scales to
many identities without slowing the common case. The first-cut libraries cover
the gamma family (`Gamma`/`LogGamma`/`PolyGamma` recurrences,
`Pochhammer`/`Beta`/`Factorial` → `Gamma`, plus guarded **pair** reductions:
reflection `Gamma[z] Gamma[1-z] -> Pi/Sin[Pi z]` for non-integer `z`, and
conjugate pairs `Gamma[1±I b] -> Pi b/Sinh[Pi b]` and
`Gamma[1/2±I b] -> Pi/Cosh[Pi b]` — so `FullSimplify[Gamma[1/4] Gamma[3/4]] ->
Pi Sqrt[2]` and `FullSimplify[Gamma[1+I] Gamma[1-I]] -> Pi Csch[Pi]`), the error
functions
(`Erf[z] + Erfc[z] -> 1`), the dilogarithm (`PolyLog[2, z] + PolyLog[2, -z] ->
PolyLog[2, z^2]/2`), real radicals (`Surd[x, n]^n -> x`), and the
**logarithm of a trig/hyperbolic cofunction combination** — the real-log
antiderivative (Gudermannian) family, which has no general algorithmic route and
so lives here rather than in `Simplify`:
`Log[Sec[x] ± Tan[x]] -> ± ArcTanh[Sin[x]]`,
`Log[Csc[x] ± Cot[x]] -> ± ArcTanh[Cos[x]]`,
`Log[Cosh[x] ± Sinh[x]] -> ± x`,
`Log[Coth[x] ± Csch[x]] -> ± ArcCoth[Cosh[x]]`, and the half-angle forms
`Log[Tan[x/2]] -> -ArcTanh[Cos[x]]`, `Log[Cot[x/2]] -> ArcTanh[Cos[x]]`,
`Log[Tanh[x/2]] -> -ArcCoth[Cosh[x]]`, `Log[Coth[x/2]] -> ArcCoth[Cosh[x]]`. Each
two-term rule is coefficient-tolerant for a provably-positive factor at either
sign (`Log[2 Sec[x] - 2 Tan[x]] -> Log[2] - ArcTanh[Sin[x]]`), matching whichever
factored or distributed form the pipeline produces; a symbolic or negative
coefficient, or an unequal pair, is left untouched (no branch is crossed).

**Options** (in addition to the positional assumption):
- `ComplexityFunction -> f` — custom complexity measure (forwarded to `Simplify`).
- `TransformationFunctions -> {f1, ...}` — extra user transforms, merged with the
  relevance-selected set (`Automatic` is always retained).
- `TimeConstraint -> {tLoc, tTot}` — at most `tLoc` seconds per individual
  transformation and `tTot` seconds in total; a bare `t` means `{t, t}` and the
  default is `Infinity`. On a total timeout the plain-`Simplify` result is
  returned, preserving the contract above.

`FullSimplify` is implemented in the Mathilda language
(`src/internal/simp/FullSimplify.m`, loaded at startup) rather than C; see
[`LoadModule`](file-io.md#loadmodule) for the runtime module-loading mechanism it
uses. `Protected`.

```mathematica
In[1]:= FullSimplify[Gamma[x + 1]/Gamma[x]]
Out[1]= x

In[2]:= FullSimplify[LogGamma[x + 1] - LogGamma[x]]
Out[2]= Log[x]

In[3]:= FullSimplify[Erf[x] + Erfc[x]]
Out[3]= 1

In[4]:= FullSimplify[PolyLog[2, z] + PolyLog[2, -z]]
Out[4]= 1/2 PolyLog[2, z^2]

In[5]:= FullSimplify[Surd[x, 3]^3]
Out[5]= x

In[6]:= FullSimplify[Pochhammer[a, n]/Gamma[a + n]]
Out[6]= 1/Gamma[a]

In[7]:= FullSimplify[Gamma[x + 1]/Gamma[x], TimeConstraint -> {1, 5}]
Out[7]= x

In[8]:= FullSimplify[Log[Sec[x] + Tan[x]]]
Out[8]= ArcTanh[Sin[x]]

In[9]:= FullSimplify[Log[Coth[x] - Csch[x]]]
Out[9]= -ArcCoth[Cosh[x]]
```

First-cut limitations: the gamma recurrence matches a literal `+1` shift only
(`Gamma[x + 2]/Gamma[x]` is left unreduced), and an identity over an `Orderless`
sum must match the whole sum (`a + Erf[x] + Erfc[x]` is not collapsed). The
identity collection is expected to grow; each family lives in its own file under
`src/internal/simp/transforms/`.

## SimplifyCount
The complexity measure used by `Simplify` when no `ComplexityFunction` option (or
`ComplexityFunction -> Automatic`) is given.

- `SimplifyCount[expr]`

**Features**:
- `Listable`, `Protected`.
- Per node: a symbol, the integer `0`, or a string counts `1`; a positive integer
  counts its decimal-digit count; a negative integer counts its digits plus one
  for the sign; `Rational[n, d]` counts `SimplifyCount[n] + SimplifyCount[d] + 1`;
  `Complex[re, im]` counts `SimplifyCount[re] + SimplifyCount[im] + 1`; a real
  (machine or MPFR) counts `2`; a function `h[a1, ...]` counts
  `SimplifyCount[h] + Sum SimplifyCount[ai]`.
- Matches Mathematica's definition, so the integer-digit penalty keeps
  `100 Log[2]` (count 6) ahead of `Log[2^100]` (count 32).

```mathematica
In[1]:= SimplifyCount[100 Log[2]]
Out[1]= 6

In[2]:= SimplifyCount[Log[2^100]]
Out[2]= 32

In[3]:= SimplifyCount[1/2]
Out[3]= 3

In[4]:= SimplifyCount[3.14]
Out[4]= 2
```

## TransformationFunctions
An option for `Simplify` giving the list of functions to apply to try to
transform parts of an expression.

- `TransformationFunctions -> Automatic` — use the built-in collection of
  transformation functions (the default).
- `TransformationFunctions -> {f1, f2, ...}` — use **only** the functions `fi`;
  the built-in pipeline is suppressed.
- `TransformationFunctions -> {Automatic, f1, ...}` — use the built-in
  transformation functions **together with** the `fi`.

**Features**:
- Each `fi` may be any function — a builtin head such as `Together` or `Cancel`,
  or a pure function such as `(# /. a -> 0 &)`.
- Every function is applied to the whole expression and to each of its
  subexpressions; the candidate of lowest complexity (per `ComplexityFunction`)
  is kept, in the same minimum-complexity search used for the built-in
  transforms.
- The option propagates through `Simplify`'s list / relation threading and
  through the inexact-input rationalise/numericalise path.

```mathematica
In[1]:= Simplify[(x^2 - 1)/(x - 1), TransformationFunctions -> {Cancel}]
Out[1]= 1 + x

In[2]:= Simplify[Sin[x]^2 + Cos[x]^2, TransformationFunctions -> {}]
Out[2]= Cos[x]^2 + Sin[x]^2
```

## Refine
Gives the form an expression would take if its symbols satisfied the assumptions.

- `Refine[expr, assum]` — gives the form of `expr` obtained if the symbols in it
  were replaced by explicit values satisfying `assum`.
- `Refine[expr]` — uses the default assumptions from any enclosing `Assuming`
  construct (`$Assumptions`).

**Options**:
- `Assumptions` (default `$Assumptions`) — default assumptions to append to
  `assum`. Given as an option value it *replaces* `$Assumptions`; the same
  positional/option policy as `Simplify` and `PossibleZeroQ`.
- `TimeConstraint` (default `30`) — seconds to spend on any single condition
  check (a `Reduce`/CAD entailment) before giving up on that transformation.

**Features**:
- `Protected`. Shares the assumption engine with `Simplify` (`AssumeCtx`,
  `apply_assumption_rules`): every rewrite `Refine` applies is one `Simplify`
  also applies. `Refine` is the rewrite-and-decide pass *without* `Simplify`'s
  complexity-minimising search, so it does not, e.g., factor.
- Assumption-driven rewrites: `Sqrt[x^2] -> x / -x / Abs[x]` (for `x > 0` /
  `x < 0` or `x <= 0` / real), `Abs[x] -> -x` for `x <= 0`,
  `(x^m)^r -> x^(m r)` (for `x >= 0`), `(a^b)^c -> a^(b c)` (for `-1 < b < 1`),
  `a^p b^p -> (a b)^p` (for `a, b > 0`), `Sign[x] -> ±1` and `Arg[x] -> 0 / Pi`
  under sign facts, `Log[x] -> I Pi + Log[-x]` (`x < 0`), `Log[x^p] -> p Log[x]`
  (`x > 0`), `Log[x^2] -> 2 Log[Abs[x]]` and `Log[E^x] -> x` (real `x`),
  `Log[x rest] -> Log[x] + Log[rest]` (positive `x`), `Sin[k Pi] -> 0` and
  `Cos[x + k Pi] -> (-1)^k Cos[x]` (integer `k`), `ArcTan[Tan[x]] -> x` on the
  principal domain, `Re`/`Im`/`Conjugate`/`Arg`/`Abs` of real-symbol expressions
  (e.g. `Abs[a + b I] -> Sqrt[a^2 + b^2]`),
  `Floor`/`Ceiling`/`Round`/`IntegerPart`/`FractionalPart` and `Mod` under
  integer / interval / modular facts. The per-symbol rule synthesis is pruned to
  the symbols the target actually uses, so a large assumption set does not
  overflow the rule buffer.
- Predicate decisions: `Element[x, dom]` via the assumption-aware domain
  prover (including compound expressions, and the sign domains
  `Positive`/`Negative`/`NonNegative`/`NonPositive`); equations via the
  assumption-aware zero test **plus** equality-substitution
  (`Refine[a == b, a - b == 0] -> True`); inequalities and their logical
  combinations via the `Reduce`/CAD entailment (`P` is `True` iff `assum && !P`
  is unsatisfiable over the reals). Quantities appearing algebraically in
  inequalities are assumed real. The zero test is sound under coupling equality
  assumptions: it never reports a genuine identity as non-zero (it downgrades to
  undecided rather than sampling points the assumptions exclude).
- Purely symbolic/structural: no packed/NDArray kernel and no `Compile[]`
  lowering (it returns symbolic expressions, not machine numbers).

```mathematica
In[1]:= Refine[Sqrt[x^2], x > 0]
Out[1]= x

In[2]:= Refine[Sqrt[x^2], Element[x, Reals]]
Out[2]= Abs[x]

In[3]:= Refine[Sign[x^2 - x y + y^2 + 1], Element[x | y, Reals]]
Out[3]= 1

In[4]:= Refine[a^2 - b^2 + 1 == 0, a + b == 0]
Out[4]= False

In[5]:= Assuming[x > 0, Refine[Sqrt[x^2 y^2], y < 0]]
Out[5]= -x y
```

## Assuming
Evaluates an expression with extra assumptions in effect.

- `Assuming[assum, expr]` — evaluates `expr` with `assum` appended to
  `$Assumptions`, so `assum` is included in the default assumptions used by
  functions such as `Simplify`.

**Features**:
- `HoldRest`, `Protected` (the assumption argument evaluates; the body is held
  until the assumption is in scope).
- Effectively `Block[{$Assumptions = $Assumptions && assum}, expr]`, so nested
  `Assuming` calls compose and the rebinding of `$Assumptions` is restored on
  exit. Lists of assumptions are converted to conjunctions.

```mathematica
In[1]:= Assuming[x > 0, Simplify[Sqrt[x^2 y^2], y < 0]]
Out[1]= -x y
```

## $Assumptions
The default setting for the `Assumptions` option used by `Simplify` and other
functions that take assumptions.

**Features**:
- A system symbol with default `OwnValue` `True` (no assumptions). `Assuming`
  temporarily extends `$Assumptions` for the duration of its body.

```mathematica
In[1]:= $Assumptions
Out[1]= True
```

## Element
Tests domain membership, consulting the current assumptions.

- `Element[x, dom]` — returns `True` if `x` is provably an element of `dom` under
  the current `$Assumptions`, `False` if it is provably not, and stays
  unevaluated otherwise.

**Features**:
- `Protected`.
- Supported domains: `Integers`, `Rationals`, `Reals`, `Algebraics`, `Complexes`,
  `Booleans`, `Primes`, `Composites`.
- Numeric and structural literals decide directly. Symbolic queries consult
  `$Assumptions`, honouring the `Integer ⊆ Rational ⊆ Algebraic ⊆ Real ⊆ Complex`
  lattice.
- `Element[{x1, ..., xN}, dom]` and `Element[x1 | ... | xN, dom]` are shorthand
  for the conjunction `Element[x1, dom] && ... && Element[xN, dom]`: `True` /
  `False` if every component decides, otherwise unevaluated (and treated as a
  joint per-variable fact by `Simplify`).

```mathematica
In[1]:= Element[7, Primes]
Out[1]= True

In[2]:= Element[5/2, Integers]
Out[2]= False

In[3]:= Element[1 + I, Reals]
Out[3]= False

In[4]:= Element[x, Reals]
Out[4]= Element[x, Reals]
```

## $SimplifyDebug
A system symbol that turns on tracing of `Simplify`'s transform pipeline.

**Features**:
- Default `False`. When set to `True`, `Simplify` prints one line per transform
  invocation to **stderr**, in the form
  `/<TransformName>/: <input> -> <output> [<elapsed> ms]`. Useful for diagnosing
  slow or hanging `Simplify` calls and runaway candidate-set growth. The value is
  read directly off the `OwnValue`, so there is no cost when it is `False`.

```mathematica
In[1]:= $SimplifyDebug = True; Simplify[a x + b x]; $SimplifyDebug = False;
(* stderr: *)
(* /PythagCanon/: a x + b x -> a x + b x [0.01 ms]                    *)
(* /TanAddition/: a x + b x -> a x + b x [0.00 ms]                    *)
(* ...                                                                *)
```
