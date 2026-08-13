# Simplify

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Simplify[expr]`**

performs a sequence of algebraic and other transformations on expr and returns the simplest form it finds.

**`Simplify[expr, assum]`**

does simplification using assumptions assum.

<details>
<summary>Notes</summary>

Options: Assumptions (default $Assumptions) -- facts assumed while simplifying. ComplexityFunction (default: leaf count plus integer-digit count, matching Mathematica) -- ranks candidate forms; the lowest-scoring form is returned. TransformationFunctions (default Automatic) -- the functions applied to try to transform parts of expr. Automatic uses the built-in collection; {f1, f2, ...} uses only the fi; {Automatic, f1, ...} uses the built-in functions together with the fi. The built-in collection tries Together, Cancel, Expand, Factor, FactorSquareFree, Apart, TrigExpand, TrigFactor, and a TrigToExp/ExpToTrig roundtrip, keeping the smallest result. Under positivity / reality assumptions Simplify also applies Log/Power identities -- Log\[a b\] -\> Log\[a\] + Log\[b\], (a b)^c -\> a^c b^c, (a^p)^q -\> a^(p q), Log\[a^p\] -\> p Log\[a\] and the like -- whenever the operand-domain conditions are provable from the assumption set. Assumptions can be equations, inequalities, domain specifications such as Element\[x, Integers\], or logical combinations of these. Lists of assumptions are converted to conjunctions. Simplify automatically threads over lists, equations, inequalities, and logic functions.

</details>

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

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

### Scope (1)

```mathematica
In[8]:= Simplify[Sqrt[2] Sqrt[3]]
Out[8]= Sqrt[6]
```

### Applications (8)

```mathematica
In[9]:= Simplify[(x^2 - 1)/(x - 1)]
Out[9]= 1 + x

In[10]:= Simplify[Sin[x]^2 + Cos[x]^2]
Out[10]= 1

In[11]:= Simplify[x + x + x]
Out[11]= 3 x

In[12]:= Simplify[Sqrt[x^2], x > 0]
Out[12]= x

In[13]:= Simplify[Sqrt[x^2], Element[x, Reals]]
Out[13]= Abs[x]

In[14]:= Simplify[Cosh[x]^2 - Sinh[x]^2]
Out[14]= 1

In[15]:= Simplify[Log[a b] - Log[a] - Log[b], {a > 0, b > 0}]
Out[15]= 0

In[16]:= Simplify[Cos[3 x]/Cos[x] - (2 Cos[2 x] - 1)]
Out[16]= 0
```

## Options & behaviour

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

A predicate that appears literally among the assumed facts folds to `True`:

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

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Discriminant of deg 20 | 2.51 s | 0.068 s | 0.182 s |
| Simplify quartic-to-Cos[4x] | 2.35 s | 0.003 s | 7.45 s |
| TrigReduce product of 4 sines | 0.894 s | 0.14 s | 16.5 s |
| Expand (1+x)^400 | 0.434 s | 0.107 s | 0.003 s |
| Cancel deg-60 over deg-58 | 0.337 s | 0.569 s | 7.37 s |
| PolynomialGCD, coprime deg 40 | 0.252 s | 0.087 s | 0.334 s |

## Implementation notes

**Algorithm.** Simplify is a complexity-weighted, memoized candidate-set search
over the existing battery of algebraic transforms. `builtin_simplify`
(`simp_builtins.c`) parses options — a positional or `Assumptions ->`
assumption (combined into the `$Assumptions` default with `And`), a
`ComplexityFunction`, and `TransformationFunctions` — builds an `AssumeCtx` of
normalised facts, then drives the search. Inexact inputs are rationalised on
entry and numericalised on exit (`internal_rationalize_then_numericalize`).
`Equal`/`Less`/`And`/... heads are threaded manually (with a relational
rebalancing candidate); `List` is threaded via `ATTR_LISTABLE`.

The core driver is `simp_bottomup` (`simp_bottomup.c`), which descends into
every `Plus`/`Times`/`Power` child (memoizing each result in a `SimpMemo` hash
table) and dispatches each node through `simp_dispatch` → `simp_search`.
`simp_classify` routes pure-polynomial/rational shapes to dedicated pipelines
(`simp_pipeline_polynomial`/`_rational`/`_logexp`) and there are top-level fast
paths: a `SHAPE_RATIONAL` shortcut that runs `Together`/`Cancel`/`Factor` once
at the top, an algebraic-tower `Together[expr, Extension -> Automatic]` collapse,
and a `simp_trig_rational` substitution that maps trig/opaque subtrees to ground
symbols and works in the quotient ring.

`simp_search` is the heuristic engine. It seeds a `CandSet` with the input plus
the output of a long list of correctness-preserving rewriters (assumption rules,
log/exp identities, `SimpLogRules`, trig roundtrip, ExpToTrig, Pythagorean
square-completion / reduction / canonicalisation, trig-at-rational-Pi, tan-
addition, half-angle, radical product combine, sqrt/cube-root/algebraic
denesting, roots-of-unity, factorial decomposition, per-variable Collect). It
then runs `SIMP_ROUNDS` (= 2) rounds in which every seed is fed through the
`SIMP_TRANSFORMS` table — `Together, Cancel, Expand, ExpandNumerator,
ExpandDenominator, Factor, FactorSquareFree, FactorTerms, Apart, TrigExpand,
TrigFactor, TrigReduce, TrigToExp` — plus chained Pythagorean/radical/trig-Pi
passes. Each transform call is gated by `transform_can_fire` (cheap precondition
check) and the running best is updated by `update_best` against the complexity
score; candidates strictly worse than their parent are dropped from the next
round (with a loosened `2*parent + 8` bound for TrigExpand and the seed-phase
blow-up guard). A short-circuit (`simp_best_is_zero`) exits as soon as a literal
`0` is reached. Final polish passes apply `simp_lift_common_factor`,
`transform_pythag_reduce`, and `canon_negate_pairs`.

**Complexity scorer.** Candidates are ranked by `score_with_func`: with no user
`ComplexityFunction` this is `simp_default_complexity` (the SimplifyCount metric
— LeafCount with integers contributing their decimal-digit count) plus a
`nested_radical_penalty` (+3 per truly-nested radical). A user function is
evaluated as `f[candidate]`; an Integer result is used directly, a BigInt scores
`SIMP_SCORE_INF`, anything else falls back to the default. Lower wins; ties
favour the form that reached the score plateau first (force-take semantics let
assumption/log/factorial rewrites win even on a tie).

**Data structures.** `CandSet` (dynamic `Expr*` array, dedup via `expr_eq`, capped
at `SIMP_CAND_CAP` = 12); `SimpMemo` (256-bucket chained hash of input→best for
the bottom-up driver); `AssumeCtx` (flat fact array from `assume_ctx_from_expr`);
a per-Simplify `FactorMemo` shared by Factor/Trig* so duplicate subexpression
work hits the cache. Each transform is invoked as `f[candidate]` through the real
evaluator (`traced_call_unary`), so Simplify composes the same builtins exposed
in the REPL.

**Limits.** Not a complete decision procedure: no real inequality reasoner
(`Simplify[x>0, x>0]` folds only by literal fact match), and assumption-driven
wins depend on the structural provers in `simp_assume.c`.

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

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Together](../../algebra/Together/), [Cancel](../../algebra/Cancel/), [Expand](../../algebra/Expand/), [Factor](../../algebra/Factor/), [FactorSquareFree](../../algebra/FactorSquareFree/), [Apart](../../algebra/Apart/), [TrigExpand](../../elementary-functions/TrigExpand/)

- Joel S. Cohen, *Computer Algebra and Symbolic Computation: Mathematical Methods* (A K Peters, 2003).
- Source: [`src/simp/simp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
- Tests: [`tests/test_assuming.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assuming.c)
- Tests: [`tests/test_cherry_dilog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog.c)
- Tests: [`tests/test_cherry_ei.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_ei.c)
- Tests: [`tests/test_cherry_li.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_li.c)

## Notes & additional examples

### Notes

`Simplify` tries a collection of transformations — `Together`, `Cancel`,
`Expand`, `Factor`, `Apart`, `TrigExpand`, `TrigFactor`, and a `TrigToExp` round
trip — and keeps the smallest result, so it can cancel `(x^2-1)/(x-1)` and reduce
the Pythagorean identity to `1`. A second argument supplies assumptions:
`Simplify[Sqrt[x^2], x > 0]` uses `x > 0` to drop the absolute value and return
`x`. Assumptions may be equations, inequalities, or domain statements like
`Element[x, Integers]`. `Simplify` threads automatically over lists, equations,
and logical combinations.
