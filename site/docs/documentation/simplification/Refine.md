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

## Examples (5)

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

- `Protected`. Shares the assumption engine with `Simplify` (`AssumeCtx`,
  `apply_assumption_rules`): every rewrite `Refine` applies is one `Simplify`
  also applies. `Refine` is the rewrite-and-decide pass *without* `Simplify`'s
  complexity-minimising search, so it does not, e.g., factor.
- Assumption-driven rewrites: `Sqrt[x^2] -> x / -x / Abs[x]`,
  `(x^m)^r -> x^(m r)` (for `x >= 0`), `(a^b)^c -> a^(b c)` (for `-1 < b < 1`),
  `a^p b^p -> (a b)^p` (for `a, b > 0`), `Log[x] -> I Pi + Log[-x]` (`x < 0`),
  `Log[x^p] -> p Log[x]` (`x > 0`), `Sin[k Pi] -> 0` and
  `Cos[x + k Pi] -> (-1)^k Cos[x]` (integer `k`), `ArcTan[Tan[x]] -> x` on the
  principal domain, `Re`/`Im`/`Conjugate` of real-symbol expressions,
  `Floor`/`Ceiling`/`Round`/`IntegerPart`/`FractionalPart` and `Mod` under
  integer / interval / modular facts.
- Predicate decisions: `Element[x, dom]` via the assumption-aware domain
  prover (including compound expressions); equations via the assumption-aware
  zero test; inequalities and their logical combinations via the `Reduce`/CAD
  entailment (`P` is `True` iff `assum && !P` is unsatisfiable over the reals).
  Quantities appearing algebraically in inequalities are assumed real.
- Purely symbolic/structural: no packed/NDArray kernel and no `Compile[]`
  lowering (it returns symbolic expressions, not machine numbers).

**Attributes:** `Protected`.

## References

**See also:** [Assuming](../../simplification/Assuming/), [$Assumptions](../../simplification/$Assumptions/), [Simplify](../../simplification/Simplify/), [PossibleZeroQ](../../expression-information/PossibleZeroQ/), [Reduce](../../solutions-of-equations/Reduce/), [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [Conjugate](../../arithmetic/Conjugate/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
- Tests: [`tests/test_refine.c`](https://github.com/stblake/mathilda/blob/main/tests/test_refine.c)
