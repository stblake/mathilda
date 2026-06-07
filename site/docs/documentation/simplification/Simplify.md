# Simplify

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Simplify[expr]
    performs a sequence of algebraic and other transformations on expr and returns the simplest form it finds.
Simplify[expr, assum]
    does simplification using assumptions assum.

Options:
  Assumptions (default $Assumptions) -- facts assumed while simplifying.
  ComplexityFunction (default: leaf count plus integer-digit count, matching Mathematica) -- ranks candidate forms; the lowest-scoring form is returned.
  TransformationFunctions (default Automatic) -- the functions applied to try to transform parts of expr. Automatic uses the built-in collection; {f1, f2, ...} uses only the fi; {Automatic, f1, ...} uses the built-in functions together with the fi.

The built-in collection tries Together, Cancel, Expand, Factor, FactorSquareFree, Apart, TrigExpand, TrigFactor, and a TrigToExp/ExpToTrig roundtrip, keeping the smallest result.
Under positivity / reality assumptions Simplify also applies Log/Power identities -- Log[a b] -> Log[a] + Log[b], (a b)^c -> a^c b^c, (a^p)^q -> a^(p q), Log[a^p] -> p Log[a] and the like -- whenever the operand-domain conditions are provable from the assumption set.
Assumptions can be equations, inequalities, domain specifications such as Element[x, Integers], or logical combinations of these. Lists of assumptions are converted to conjunctions.
Simplify automatically threads over lists, equations, inequalities, and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

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

```mathematica
In[1]:= Simplify[Sqrt[2] Sqrt[3]]
Out[1]= Sqrt[6]
```

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

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Joel S. Cohen, *Computer Algebra and Symbolic Computation: Mathematical Methods* (A K Peters, 2003).
- Source: [`src/simp/simp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Simplify[(x^2 - 1)/(x - 1)]
Out[1]= 1 + x
```

```mathematica
In[1]:= Simplify[Sin[x]^2 + Cos[x]^2]
Out[1]= 1
```

```mathematica
In[1]:= Simplify[x + x + x]
Out[1]= 3 x
```

```mathematica
In[1]:= Simplify[Sqrt[x^2], x > 0]
Out[1]= x
```

### Notes

`Simplify` tries a collection of transformations — `Together`, `Cancel`,
`Expand`, `Factor`, `Apart`, `TrigExpand`, `TrigFactor`, and a `TrigToExp` round
trip — and keeps the smallest result, so it can cancel `(x^2-1)/(x-1)` and reduce
the Pythagorean identity to `1`. A second argument supplies assumptions:
`Simplify[Sqrt[x^2], x > 0]` uses `x > 0` to drop the absolute value and return
`x`. Assumptions may be equations, inequalities, or domain statements like
`Element[x, Integers]`. `Simplify` threads automatically over lists, equations,
and logical combinations.
