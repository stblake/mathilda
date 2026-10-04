# Refine

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Refine[expr, assum]`**

gives the form of expr that would be obtained if the symbols in it were replaced by explicit values satisfying the assumptions assum.

**`Refine[expr]`**

uses the default assumptions specified by any enclosing Assuming constructs ($Assumptions).

<details>
<summary>Notes</summary>

Options: Assumptions (default $Assumptions) -- default assumptions to append to assum. Given as an option value, it prevents Refine from also using $Assumptions. TimeConstraint (default 30) -- how many seconds to spend on any single condition check before giving up on that transformation. Assumptions can be equations, inequalities, domain specifications such as Element\[x, Integers\], or logical combinations of these. Quantities appearing algebraically in inequalities are assumed to be real. Refine can be applied to expressions, and to equations, inequalities, and domain specifications, which it decides to True or False when they provably follow from (or contradict) the assumptions. Refine is one of the transformations tried by Simplify; use Simplify or FullSimplify for a broader search.

</details>

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

### Applications (9)

The square root of a square collapses once the sign is known

```mathematica
In[6]:= Refine[Sqrt[x^2], x > 0]
Out[6]= x
```

Abs picks up the sign fact

```mathematica
In[7]:= Refine[Abs[x], x < 0]
Out[7]= -x
```

Log[x^p] -> p Log[x] for positive x

```mathematica
In[8]:= Refine[Log[x^2], x > 0]
Out[8]= 2 Log[x]
```

Sign resolved to +1 under the assumption

```mathematica
In[9]:= Refine[Sign[x], x > 0]
Out[9]= 1
```

A real symbol is its own conjugate

```mathematica
In[10]:= Refine[Conjugate[x], Element[x, Reals]]
Out[10]= x
```

An inequality decided by a Reduce/CAD entailment

```mathematica
In[11]:= Refine[x > 0, x > 1]
Out[11]= True
```

Algebraic quantities in an inequality are assumed real

```mathematica
In[12]:= Refine[Element[x, Reals], x > 0]
Out[12]= True
```

Equality substitution drives the difference to zero

```mathematica
In[13]:= Refine[a - b, a == b]
Out[13]= 0
```

Refine with no positional assumption reads $Assumptions

```mathematica
In[14]:= Assuming[a > 0 && b > 0, Refine[Sqrt[a^2 b^2]]]
Out[14]= a b
```

## Algorithm

Mathilda -- Refine implementation.

Refine[expr, assum] gives the form of `expr` that would be obtained if the symbols in it were replaced by explicit values satisfying `assum`. It is a thin orchestrator over machinery that already exists and was deliberately shared for it (simp.h reserves AssumeCtx "so future modules (Refine, ...) can share it"):

```text
  1. Argument / option parsing mirrors builtin_possible_zero_q
     (src/zero_test.c): one positional assumption, an Assumptions -> X option
     that overrides $Assumptions, and a TimeConstraint option. The effective
     assumption is And-combined exactly as PossibleZeroQ/Simplify do.

  2. Predicate expressions are decided to True/False:
       Element[x, dom]              -> element_decide (assumption-aware,
                                       including compound-expression domain
                                       inference via prov_re/prov_int).
       Equal / Unequal              -> zero_test_decide_assuming on lhs - rhs,
                                       with a Reduce fallback.
       Less / ... / logic           -> Reduce/CAD entailment over the reals:
                                       P is True  iff Reduce[A && !P] is False,
                                       P is False iff Reduce[A &&  P] is False.

  3. Every other expression is rewritten by apply_assumption_rules
     (Sqrt[x^2]->x/-x/Abs[x], Log[x^p]->p Log[x], integer-k trig, Abs/Sign/
     Conjugate under sign facts, Floor/Ceiling/Mod under integer/interval
     facts, ...) followed by a deep-positivity post-pass that resolves
     Sign[p]/Abs[p]/Sqrt[p^2] for symbolic polynomials p that the fast
     sign prover cannot settle, using the same Reduce entailment.
```

TimeConstraint (default 30s) is enforced with Simplify's cooperative wall-clock budget (simp_set_time_budget / simp_mono_seconds): Reduce/CAD has no interior abort hook, so the deadline is checked BEFORE each entailment call and a running Reduce is never preempted (best-effort, and never the malloc-lock-prone async TimeConstrained[]).

Refine is a purely symbolic/structural head: it has no element-wise numeric semantics, so it deliberately carries no NDArray/packed kernel and no Compile[] lowering.

Ownership follows the builtin contract: return a new tree or steal from

```text
`res`; never expr_free(res).
```

## Implementation notes

**Algorithm.** `builtin_refine` is a thin orchestrator over the assumption engine
`Simplify` already uses (`AssumeCtx`, `apply_assumption_rules`, declared in `simp.h`,
which reserves that type "so future modules (Refine, …) can share it"). It first splits
trailing options from positional arguments and forms the *effective* assumption with the
identical policy to `PossibleZeroQ`/`Simplify`: an `Assumptions -> X` option replaces
`$Assumptions`, a positional assumption is `And`-conjoined with `$Assumptions`
(`read_dollar_assumptions`), and the result is evaluated so `And[True, p]` canonicalises to
`p`. That expression is parsed into an `AssumeCtx` fact set by `assume_ctx_from_expr`. From
there two paths run. **(1) Predicate head** (`Element`, or a relational/logical head —
`Equal`, `Less`, `And`, …): `refine_decide_predicate` tries to settle it to `True`/`False`.
`Element[x, dom]` goes to the assumption-aware `element_decide`; `Equal`/`Unequal` first
run `zero_test_decide_assuming` on `lhs − rhs` (and an equality-substitution via
`apply_assumption_rules`), then fall back to `Reduce` over the complexes; inequalities and
logic go to a `Reduce`/CAD entailment over the reals, where `P` is `True` iff `Reduce[A &&
!P]` is unsatisfiable (literally `False`) and `False` iff `Reduce[A && P]` is. **(2) Rewrite
path** (everything else): `apply_assumption_rules` rewrites the expression under the facts
(the very rules `Simplify` applies — `Sqrt[x^2]->±x/Abs[x]`, `Log[x^p]->p Log[x]`, integer-`k`
trig, `Abs`/`Sign`/`Conjugate` under sign facts, …), then a bottom-up `deep_positivity_walk`
resolves `Sign[p]`, `Abs[p]` and `Sqrt[p^2]` for compound polynomials `p` the fast sign
prover could not settle, each via a `deep_sign` CAD run, and the tree is evaluated. With no
usable facts the call is the identity (it steals the positional expression out of `res`).

**Data structures.** Everything is `Expr` trees. Assumptions are held as an `AssumeCtx`
(flat `Expr*` fact array) borrowed from `simp`. Entailment is threaded through a
`RefineBudget { double deadline; int calls_left; }`; each `Reduce` query is built as
`Reduce[stmt, varlist, domain]`, where `collect_bare_vars` gathers the distinct bare-symbol
leaves (skipping protected real constants such as `Pi`, and *not* reusing the polynomial
`collect_variables`, which would atomise `x^p`/`Floor[x]` into pseudo-variables `Reduce`
rejects).

**Complexity / limits.** Dominated by the `Reduce`/CAD entailment, which is doubly
exponential in the number of variables; it is bounded two ways — the variable set is capped
at 6 (wider statements decline) and at most `REFINE_MAX_ENTAILMENT_CALLS` (24) `Reduce`
calls run per `Refine`. `TimeConstraint` (default 30 s) is a cooperative wall-clock budget
checked *before* each entailment via `simp_mono_seconds`; a running `Reduce` is never
preempted (CAD has no interior abort hook), and the async, malloc-lock-prone
`TimeConstrained[]` is deliberately avoided. Predicates it cannot decide, and expressions no
rule rewrites, pass through unchanged. `Refine` is purely symbolic/structural: it carries no
NDArray/packed kernel and no `Compile[]` lowering, by design, since it returns symbolic
expressions rather than machine numbers.

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

**Attributes:** `Protected`.

## References

**See also:** [Assuming](../../simplification/Assuming/), [$Assumptions](../../simplification/$Assumptions/), [Simplify](../../simplification/Simplify/), [PossibleZeroQ](../../expression-information/PossibleZeroQ/), [Reduce](../../solutions-of-equations/Reduce/), [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [Conjugate](../../arithmetic/Conjugate/)

- Source: [`src/refine.c`](https://github.com/stblake/mathilda/blob/main/src/refine.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
- Tests: [`tests/test_refine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_refine.c)

## Notes & additional examples

### Notes

`Refine[expr, assum]` gives the form `expr` would take if its symbols were replaced by
explicit values satisfying `assum`. It shares its whole assumption engine with `Simplify`
(`AssumeCtx`, `apply_assumption_rules`): every rewrite `Refine` applies is one `Simplify`
also applies. The difference is that `Refine` is the rewrite-and-decide pass *without*
`Simplify`'s complexity-minimising search — so it will turn `Sqrt[x^2]` into `x` under
`x > 0`, but it does not, for instance, factor.

A predicate target (`Element[...]`, an equation, an inequality, or a logical combination) is
decided to `True`/`False` when it provably follows from or contradicts the assumptions:
domain membership through the assumption-aware prover, equations through the assumption-aware
zero test plus equality substitution, and inequalities through a `Reduce`/CAD entailment over
the reals (`P` is `True` exactly when `assum && !P` is unsatisfiable). Quantities appearing
algebraically in inequalities are assumed real.

`Refine[expr]` with no positional assumption uses the default assumptions — the current
`$Assumptions`, as extended by any enclosing `Assuming`. An explicit `Assumptions -> X`
option replaces `$Assumptions`; a positional assumption is conjoined with it. The
`TimeConstraint` option (default 30 s) caps the time spent on any single condition check
before that transformation is abandoned.
