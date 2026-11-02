# TrigReduce Implementation Plan

Status: planning. No code written yet. This document is the design for a
new builtin `TrigReduce[expr]` and the matching `simp_search` integration
that lets `Simplify` recognise angle-addition forms as the canonical
"reduced" target.

## 1. What TrigReduce does

`TrigReduce` is the inverse of `TrigExpand`. Where `TrigExpand` rewrites
sums and integer multiples in trig arguments into products of single-
argument trig calls, `TrigReduce` rewrites *products and powers* of
single-argument trig calls into single trig calls of *sums* of arguments.

It is the standard product-to-sum / power-reduction pass:

| Input                         | TrigReduce output                          |
| ----------------------------- | ------------------------------------------ |
| `Sin[a] Cos[b]`               | `(Sin[a + b] + Sin[a - b]) / 2`            |
| `Cos[a] Cos[b]`               | `(Cos[a + b] + Cos[a - b]) / 2`            |
| `Sin[a] Sin[b]`               | `(Cos[a - b] - Cos[a + b]) / 2`            |
| `Sin[x]^2`                    | `(1 - Cos[2 x]) / 2`                       |
| `Cos[x]^2`                    | `(1 + Cos[2 x]) / 2`                       |
| `Sin[x]^3`                    | `(3 Sin[x] - Sin[3 x]) / 4`                |
| `Cos[x]^3`                    | `(3 Cos[x] + Cos[3 x]) / 4`                |
| `Sin[a] Cos[b] + Cos[a] Sin[b]` | `Sin[a + b]`                             |
| `Cos[a] Cos[b] - Sin[a] Sin[b]` | `Cos[a + b]`                             |

Hyperbolic analogues:

| Input                          | TrigReduce output                          |
| ------------------------------ | ------------------------------------------ |
| `Sinh[a] Cosh[b]`              | `(Sinh[a + b] + Sinh[a - b]) / 2`          |
| `Cosh[a] Cosh[b]`              | `(Cosh[a + b] + Cosh[a - b]) / 2`          |
| `Sinh[a] Sinh[b]`              | `(Cosh[a + b] - Cosh[a - b]) / 2`          |
| `Sinh[x]^2`                    | `(-1 + Cosh[2 x]) / 2`                     |
| `Cosh[x]^2`                    | `(1 + Cosh[2 x]) / 2`                      |
| `Sinh[a] Cosh[b] + Cosh[a] Sinh[b]` | `Sinh[a + b]`                         |
| `Cosh[a] Cosh[b] + Sinh[a] Sinh[b]` | `Cosh[a + b]`                         |

The result is a polynomial in single-trig-of-sum-of-arguments terms,
with integer or rational coefficients.

## 2. Motivation: failing Simplify cases

Currently Mathilda has `TrigExpand` (forward direction) but no inverse.
Mathematica's `Simplify` uses `TrigReduce` internally to recognise
angle-addition forms as more compact than their expanded equivalents.
Concrete failures from the user's recent test list:

- **#20**: `Simplify[Sin[a] (Cos[b] - Sin[b]) + Cos[a] (Sin[b] + Cos[b])]`
  - Currently: returned unchanged.
  - Mathematica: `Cos[a + b] + Sin[a + b]`.
  - The expanded form has 12 leaves; the reduced form has 8. Once
    `TrigReduce` exists, `simp_search` will see the reduced form as a
    candidate and the leaf-count tiebreak picks it.
- **#13**: `Simplify[Sin[x] / Sin[x/2]]` should yield `2 Cos[x/2]`.
  - The path: `TrigReduce[Sin[x]] -> 2 Sin[x/2] Cos[x/2]` is the wrong
    direction (TrigReduce takes products to sums, not the other way). A
    proper fix needs an *inverse* power-reduction step too: half-angle
    expansion of a single trig at a doubled argument. That belongs in a
    separate `TrigReduce` mode or in `TrigExpand`'s multiple-angle path
    rather than here -- listed as out-of-scope below.
- **#12**: `Simplify[Cos[x/6] / Sin[x/3]]` should yield `1/2 Csc[x/6]`.
  - Same out-of-scope reason as #13: requires reducing `Sin[x/3]` to
    `2 Sin[x/6] Cos[x/6]`, which is a *forward* direction.
- **#11**: `Simplify[Sin[x/3] Sin[2 x] Sin[x/4] / (Cos[x/6] Cos[x/8])]`
  should yield `8 Cos[x] Sin[x/8] Sin[x/6] Sin[x]`.
  - This combines TrigReduce (collapse Sin Sin pairs to Cos differences)
    with multiple-angle expansion (Sin[2 x] -> 2 Sin[x] Cos[x]). Partial
    win: TrigReduce alone doesn't get all the way to Mathematica's form,
    but it gets significantly closer than today.

So **TrigReduce primarily fixes #20**, opens the door to more aggressive
simplification of trig-product sums, and is a prerequisite for any
bidirectional trig identity work. #11-#13 also need orthogonal
double/half-angle work that is out of scope for this plan.

## 3. Algorithm

Rule-based, in the same spirit as `TrigExpand` (which is a
`ReplaceRepeated` over a parsed rule list inside `trigsimp.c`).

### 3.1 Core rule list

Five product-to-sum rules (and their hyperbolic analogues), plus the
even-power reductions, plus the angle-addition recognition rules:

```mathematica
{
  (* Circular product-to-sum *)
  Sin[a_] Cos[b_]                    :> (Sin[a + b] + Sin[a - b]) / 2,
  Cos[a_] Sin[b_]                    :> (Sin[a + b] - Sin[a - b]) / 2,
  Sin[a_] Sin[b_]                    :> (Cos[a - b] - Cos[a + b]) / 2,
  Cos[a_] Cos[b_]                    :> (Cos[a + b] + Cos[a - b]) / 2,

  (* Hyperbolic product-to-sum *)
  Sinh[a_] Cosh[b_]                  :> (Sinh[a + b] + Sinh[a - b]) / 2,
  Cosh[a_] Sinh[b_]                  :> (Sinh[a + b] - Sinh[a - b]) / 2,
  Sinh[a_] Sinh[b_]                  :> (Cosh[a + b] - Cosh[a - b]) / 2,
  Cosh[a_] Cosh[b_]                  :> (Cosh[a + b] + Cosh[a - b]) / 2,

  (* Power reduction (degree 2) *)
  Sin[x_]^2                          :> (1 - Cos[2 x]) / 2,
  Cos[x_]^2                          :> (1 + Cos[2 x]) / 2,
  Sinh[x_]^2                         :> (-1 + Cosh[2 x]) / 2,
  Cosh[x_]^2                         :> (1 + Cosh[2 x]) / 2,

  (* Higher even/odd powers handled by repeated descent:
     Sin[x]^n gets one Sin[x]^2 pulled off, becomes
     Sin[x]^(n-2) (1 - Cos[2 x]) / 2, and the rule fires again on the
     leftover Sin[x]^(n-2). Achieved automatically by ReplaceRepeated
     because the n=2 rule is the base case. *)

  (* Angle-addition recognition.
     These are the IDEMPOTENT collapse rules: applying the four product-
     to-sum rules above to (Sin[a+b]+Sin[a-b])/2 + (Sin[a+b]-Sin[a-b])/2
     -- i.e. Sin[a] Cos[b] + Cos[a] Sin[b] -- recovers Sin[a+b]; we run
     them as a fixed-point cleanup AFTER product-to-sum so the user-
     facing TrigReduce[Sin[a] Cos[b] + Cos[a] Sin[b]] returns Sin[a+b]
     directly without leaving 2*Sin[a+b]/2 littered around. The rule
     list is short enough to encode directly. *)
  Sin[a_] Cos[b_] + Cos[a_] Sin[b_]  :> Sin[a + b],
  Sin[a_] Cos[b_] - Cos[a_] Sin[b_]  :> Sin[a - b],
  Cos[a_] Cos[b_] - Sin[a_] Sin[b_]  :> Cos[a + b],
  Cos[a_] Cos[b_] + Sin[a_] Sin[b_]  :> Cos[a - b],
  Sinh[a_] Cosh[b_] + Cosh[a_] Sinh[b_] :> Sinh[a + b],
  Sinh[a_] Cosh[b_] - Cosh[a_] Sinh[b_] :> Sinh[a - b],
  Cosh[a_] Cosh[b_] + Sinh[a_] Sinh[b_] :> Cosh[a + b],
  Cosh[a_] Cosh[b_] - Sinh[a_] Sinh[b_] :> Cosh[a - b]
}
```

### 3.2 Pipeline

```c
Expr* builtin_trigreduce_impl(Expr* res) {
    arg = res->data.function.args[0];

    // Threading: List + the same set TrigExpand threads over
    // (Equal/Less/And/Or/...). Reuse trigexpand_threads_over.

    // Step 1: ReplaceRepeated against the product-to-sum rules.
    // This brings every product/power of single-arg trig into a sum of
    // single trig calls of compound arguments.
    intermediate = ReplaceRepeated[arg, trig_reduce_rules];

    // Step 2: Expand. After product-to-sum, the result is typically a
    // rational with a constant denominator (powers of 2). Expand
    // distributes the 1/2^k constants over the sum so like-argument
    // trig calls coalesce.
    intermediate = Expand[intermediate];

    // Step 3: ReplaceRepeated against the angle-addition collapse
    // rules. Idempotent on inputs that don't have the matching shape.
    result = ReplaceRepeated[intermediate, trig_reduce_collapse_rules];

    // Step 4: Together to renormalise the constant denominators that
    // didn't get cleared by Expand.
    return Together[result];
}
```

### 3.3 Memoisation

Reuse the existing `trig_memo_call` wrapper (already used for
`TrigExpand` and `TrigFactor`). Key under the head `"TrigReduce"` so it
doesn't collide with sibling memo entries.

## 4. Module placement

All in `src/trigsimp.c`. Follows the exact pattern of `TrigExpand`:

- Two static rule trees: `trig_reduce_rules`, `trig_reduce_collapse_rules`.
  Built once in `trigsimp_init` from string literals using
  `parse_expression`.
- `builtin_trigreduce_impl(Expr*)` implements the four-step pipeline.
- `builtin_trigreduce(Expr*)` wraps it in `trig_memo_call`.
- `trigsimp_init` registers the builtin with attributes
  `ATTR_LISTABLE | ATTR_PROTECTED` and a docstring (description +
  examples).

## 5. Integration with simp_search

Two integration points in `src/simp.c`:

### 5.1 Add to `SIMP_TRANSFORMS`

Add `"TrigReduce"` as a new entry in the `SIMP_TRANSFORMS[]` list at
~line 2326 in `simp.c`. The round-loop dispatcher (`simp_search`) will
then try `TrigReduce` on every seed and propagate the result if the
score is no greater than the parent's.

The existing `transform_can_fire("TrigReduce", e, ctx)` predicate
should require:
- The input contains at least one trig or hyperbolic head
  (reuse `contains_trig_or_hyperbolic`).
- The input contains a Plus or Times somewhere (otherwise there's
  no product-to-sum work to do; TrigReduce on a single trig call is a
  no-op).

This second guard is important: without it, every Sin/Cos/Sinh/Cosh
seed in the search would fire TrigReduce and waste a memo slot per
single-trig leaf.

### 5.2 Score-gate behaviour

Default score-gate behaviour is appropriate. Unlike `TrigExpand`
(which typically grows the leaf count), `TrigReduce` typically
*shrinks* it when the input has the angle-addition shape (#20: 12 -> 8
leaves) and is roughly leaf-count-neutral otherwise. So no special
bypass is needed; the standard round-loop gate (propagate iff
`score(r) <= parent_score`) handles it correctly.

### 5.3 Chained TrigReduce on transform outputs?

Optional. The chained `simp_radicals` and `transform_pythag_reduce`
passes in the round loop exist because their inputs are produced by
`Together` / `Cancel` / etc. and would otherwise not be re-seeded in
the final round. `TrigReduce`'s analogous case would be: a `Together`
output that surfaces a `Sin[a] Cos[b] + Cos[a] Sin[b]` shape that
collapses to `Sin[a + b]`. **Unclear whether this is worth the per-
candidate cost** -- defer until we see a real test case that needs it,
to avoid speculative complexity.

## 6. Edge cases

### 6.1 Idempotent on already-reduced inputs

`TrigReduce[Sin[a + b]]` should return `Sin[a + b]` unchanged. The
rules don't fire on a single trig of a sum, only on products of
single-argument trigs, so the input passes through both
`ReplaceRepeated`s untouched. Verified by the structure of the rules.

### 6.2 Inputs with no trig

`TrigReduce[x^2 + 1]` should return `x^2 + 1` unchanged. The
ReplaceRepeated finds no match; Expand and Together leave it alone.
Cheap.

### 6.3 Mixed circular / hyperbolic

`TrigReduce[Sin[x] Cosh[y]]` has no product-to-sum identity (the
identities only apply within the same family). The rules don't fire,
the result is the input. Correct.

### 6.4 Trig of zero

`TrigReduce[Sin[a] Cos[a]]` should yield `Sin[2 a] / 2`, not `Sin[2 a] / 2 + Sin[0] / 2`.
The auto-evaluator handles `Sin[0] -> 0` and `Cos[0] -> 1` already, so
the `(Sin[a + a] + Sin[a - a]) / 2 = Sin[2 a] / 2 + 0 / 2` reduces
inside `Expand`. Correct.

### 6.5 Threading

Match `TrigExpand`'s `trigexpand_threads_over` set: List, Equal,
Unequal, Less, LessEqual, Greater, GreaterEqual, SameQ, UnsameQ, And,
Or, Not, Xor, Implies. List is handled via `ATTR_LISTABLE`; the rest
need explicit threading inside `builtin_trigreduce_impl`. Refactor
`trigexpand_threads_over` -> generic `trig_threads_over` and reuse.

### 6.6 Negative integer power

`TrigReduce[Sin[x]^(-2)]` -- not a target. The product-to-sum identity
only applies to positive integer powers. The rules use `Sin[x_]^2`
literally; the `n=-2` case won't match. Correct.

### 6.7 Symbolic exponent

`TrigReduce[Sin[x]^n]` with symbolic `n` -- no rule fires, returned
unchanged. Correct (we cannot reduce a symbolic-power trig).

## 7. Out of scope (not in this plan)

- **Half-angle expansion**: `Sin[x] -> 2 Sin[x/2] Cos[x/2]` (the
  reverse direction). This is what cases #11-#13 actually need.
  Belongs in `TrigExpand`'s multiple-angle path or a separate
  `HalfAngle` mode -- distinct work item.
- **Tan/Cot product-to-sum**: Mathematica's `TrigReduce` rewrites
  `Tan` and `Cot` via `Sin`/`Cos`. Implementable as a pre-pass that
  rewrites Tan/Cot to Sin/Cos before the main rules, then folds back.
  Skipped for the first cut to keep the rule list tight; can be added
  as a follow-on once the core is solid.
- **Chained TrigReduce on round-loop transform outputs**: see §5.3.
- **TrigReduce with explicit variable list**: Mathematica supports
  `TrigReduce[expr, vars]` to restrict reduction to specific
  arguments. Not needed for Simplify integration; can be added later
  if user-facing requests come in.

## 8. Tests

Add to `tests/test_trigsimp.c` (creating it if it doesn't exist; the
existing trig tests live in `test_trig.c` and `test_trigfactor.c`).

### 8.1 Direct TrigReduce tests

```c
test_trigreduce_sin_cos_product
    TrigReduce[Sin[a] Cos[b]] == (Sin[a + b] + Sin[a - b]) / 2

test_trigreduce_cos_cos_product
    TrigReduce[Cos[a] Cos[b]] == (Cos[a + b] + Cos[a - b]) / 2

test_trigreduce_sin_sin_product
    TrigReduce[Sin[a] Sin[b]] == (Cos[a - b] - Cos[a + b]) / 2

test_trigreduce_sin_squared
    TrigReduce[Sin[x]^2] == (1 - Cos[2 x]) / 2

test_trigreduce_cos_cubed
    TrigReduce[Cos[x]^3] == (3 Cos[x] + Cos[3 x]) / 4

test_trigreduce_pythagorean_collapse
    TrigReduce[Sin[x]^2 + Cos[x]^2] == 1

test_trigreduce_double_angle_recognition
    TrigReduce[2 Sin[x] Cos[x]] == Sin[2 x]

test_trigreduce_angle_addition_sin
    TrigReduce[Sin[a] Cos[b] + Cos[a] Sin[b]] == Sin[a + b]

test_trigreduce_angle_addition_cos
    TrigReduce[Cos[a] Cos[b] - Sin[a] Sin[b]] == Cos[a + b]

test_trigreduce_hyperbolic_angle_addition
    TrigReduce[Sinh[a] Cosh[b] + Cosh[a] Sinh[b]] == Sinh[a + b]

test_trigreduce_no_op_on_atoms
    TrigReduce[Sin[x]] == Sin[x]
    TrigReduce[x^2 + 1] == x^2 + 1
    TrigReduce[Sin[a + b]] == Sin[a + b]

test_trigreduce_threads_over_list
    TrigReduce[{Sin[x]^2, Cos[x]^2}] == {(1 - Cos[2 x]) / 2,
                                          (1 + Cos[2 x]) / 2}

test_trigreduce_threads_over_equation
    TrigReduce[Sin[x]^2 + Cos[x]^2 == y] == (1 == y)
```

### 8.2 Simplify-via-TrigReduce regression tests

In `tests/test_simplify.c`:

```c
test_simplify_angle_addition_sin_cos_sum   // covers #20
    Simplify[Sin[a] (Cos[b] - Sin[b]) + Cos[a] (Sin[b] + Cos[b])]
        == Cos[a + b] + Sin[a + b]

test_simplify_double_angle_via_trigreduce
    Simplify[2 Sin[x] Cos[x]] == Sin[2 x]      // already passes;
                                                // regression guard
```

### 8.3 Performance check

Add a Timing-based assertion (or just a manual smoke test in the
PR description, since timing is host-dependent) that
`Simplify[Sin[a] (Cos[b] - Sin[b]) + Cos[a] (Sin[b] + Cos[b])]`
finishes in under 1 second on a typical workstation.

## 9. Documentation

- Update `Mathilda_spec.md` with a `TrigReduce` entry under the
  trigonometric-simplification section, plus a "Simplify integration"
  paragraph noting the new transform in the round-loop dispatcher.
- Set the docstring via `symtab_set_docstring("TrigReduce", ...)` with
  a short description and the four core examples (Sin*Cos, Sin^2,
  angle-addition, threading over List).

## 10. Implementation steps

1. Refactor `trigexpand_threads_over` -> `trig_threads_over` (rename +
   move to a place where both functions can use it).
2. Add the two static rule trees and the parse_expression
   initialisation in `trigsimp_init`.
3. Implement `builtin_trigreduce_impl` and `builtin_trigreduce`,
   mirroring `builtin_trigexpand`'s structure. Wire memoisation.
4. Register `TrigReduce` in `trigsimp_init` with attributes and
   docstring.
5. Add `"TrigReduce"` to `SIMP_TRANSFORMS[]` in `simp.c`. Add the
   `transform_can_fire` guard for trig-presence + Plus-or-Times.
6. Write the direct tests in `tests/test_trigsimp.c` (or extend
   `test_trig.c`). Wire into `tests/CMakeLists.txt`.
7. Add the Simplify regression test for case #20 in
   `tests/test_simplify.c`.
8. Run the full test suite; verify no regressions vs. the
   pre-existing 14 FAIL lines in 6 files
   (core, integrals, linalg, parfrac, polymod, trigfactor).
9. Update `Mathilda_spec.md`.
10. Commit and push.

Estimated total: ~150-200 lines of C plus the rule strings and tests.
Most of the volume is in mirroring the existing `TrigExpand` plumbing;
the new logic is contained in the two rule lists and the four-step
pipeline body.
