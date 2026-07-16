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
- **Radical denesting** — `Sqrt[A + Sqrt[B]]` and cube-root towers collapse via
  the half-sum identity when the result is cleaner.
- **Inverse trig / hyperbolic identities** — standard relations such as
  `Sin[ArcCos[x]] == Sqrt[1 - x^2]` and `ArcSin[x] + ArcCos[x] == Pi/2` reduce.
- **Logarithm simplification** — `Log` of a positive rational is decomposed over
  its prime factors, and linear combinations of logs are fused
  (`Sum c_i Log[a_i] -> Log[Prod a_i^c_i]`).
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
  `Sin[x] Cos[y] + Cos[x] Sin[y] - Sin[x + y] -> 0` collapse too.
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
  or more distinct *positive* radical bases (e.g. `a^(1/3)` and `(a+b x)^(1/3)`)
  is reduced in the quotient ring `K[g_1, ..., g_n] / <g_k^{q_k} - base_k>`: each
  base is carried as an algebraic generator, the terms are combined over a common
  denominator, reduced modulo the generator relations, and the denominator is
  rationalised, before the radicals are substituted back. This recovers
  cross-base cancellations the single-generator `Together`/`Cancel` path cannot —
  e.g. `D[Integrate[1/(x^3 (a+b x)^(1/3)), x], x] // Simplify ->
  1/(x^3 (a+b x)^(1/3))`. Bases that are not provably positive (negative or
  complex numeric radicands) are left untouched for branch-cut safety, and the
  result is adopted only when its `SimplifyCount` strictly improves.
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
PolyLog[2, z^2]/2`), and real radicals (`Surd[x, n]^n -> x`).

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
