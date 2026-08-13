# TrigFactor

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TrigFactor[expr]`**

factors trigonometric functions in expr. TrigFactor operates on both circular and hyperbolic functions. TrigFactor factors polynomials in trigonometric functions and collapses Pythagorean, angle-addition, and double-angle identities where possible, broadly acting as the inverse of TrigExpand. TrigFactor automatically threads over lists, as well as equations, inequalities, and logic functions.

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (8)

```mathematica
In[1]:= TrigFactor[Sin[x]^2 + Cos[x]^2]
Out[1]= 1

In[2]:= TrigFactor[Cosh[x]^2 - Sinh[x]^2]
Out[2]= 1

In[3]:= TrigFactor[2 Sin[x] Cos[x]]
Out[3]= Sin[2 x]

In[4]:= TrigFactor[Cos[x]^2 - Sin[x]^2]
Out[4]= Cos[2 x]

In[5]:= TrigFactor[Sin[a] Cos[b] + Cos[a] Sin[b]]
Out[5]= Sin[a + b]

In[6]:= TrigFactor[Cos[a] Cos[b] + Sin[a] Sin[b]]
Out[6]= Cos[a - b]

In[7]:= TrigFactor[Sin[x]^2 + Tan[x]^2]
Out[7]= (1 + Cos[x]^2) Tan[x]^2

In[8]:= TrigFactor[Cosh[x]^2 - Cosh[x]^4]
Out[8]= -Cosh[x]^2 Sinh[x]^2
```

### Applications (5)

```mathematica
In[9]:= TrigFactor[Sin[a] Cos[b] + Cos[a] Sin[b]]
Out[9]= Sin[a + b]

In[10]:= TrigFactor[Cos[a] Cos[b] - Sin[a] Sin[b]]
Out[10]= Cos[a + b]

In[11]:= TrigFactor[Sin[x]^2 - Cos[x]^2]
Out[11]= -Cos[2 x]

In[12]:= TrigFactor[Sin[x]^2 + 2 Sin[x] Cos[x] + Cos[x]^2]
Out[12]= 2 Sin[1/4 Pi + x]^2

In[13]:= TrigFactor[Sinh[x]^2 + Cosh[x]^2]
Out[13]= Cosh[2 x]
```

## Implementation notes

**Algorithm.** `builtin_trigfactor_impl` factors trig expressions — broadly the
inverse of TrigExpand for the structural identities both support. It tries two
parallel paths and keeps the shorter result (by `trigfactor_leaf_count`):

- **Path A:** `trigfactor_run_pipeline` on the argument as-is. If this changes
  the expression at all it is trusted and Path B is skipped (Path B can be
  expensive on angle-sum arguments).
- **Path B:** only attempted when Path A was a no-op *and* the input has compound
  trig structure (`has_compound_trig_structure` — a `Power[trig, k≥2]` or a
  `Times` of two or more trig atoms). First `TrigExpand` the argument, then run
  the pipeline; kept only if strictly shorter than Path A. This catches
  cancellations that surface only after angle-sum expansion (e.g. `Cos[x+y] +
  Sin[x] Sin[y] -> Cos[x] Cos[y]`).

`trigfactor_run_pipeline` (trig canonicalizer suppressed throughout):
1. `ReplaceRepeated` with `trig_factor_to_sincos` — rewrite reciprocal heads as
   `Sin/Cos` ratios so `Factor` sees full polynomial structure.
2. `Together` over a common denominator.
3. `Factor` — Mathilda `Factor` treats trig atoms as polynomial variables.
   Skipped (left as the `Together` output) when the polynomial would stall:
   more than `TRIG_FACTOR_ATOM_THRESHOLD` distinct squared atoms, max trig-atom
   power above `TRIG_FACTOR_DEGREE_THRESHOLD`, or more than
   `TRIG_FACTOR_TOTAL_ATOM_THRESHOLD` distinct atoms.
4. `ReplaceRepeated` with `trig_factor_identities` — Pythagorean collapses (both
   signs, circular and hyperbolic), reverse angle-addition, reverse double-angle,
   factored-form `(Cos-Sin)(Cos+Sin) -> Cos[2x]` and `(Cosh±Sinh)` collapses, and
   the Weierstrass-style linear-combination factoring `a Sin[x] + b Cos[x] ->
   Sqrt[a^2+b^2] Sin[x + ArcTan[a, b]]` (gated by `NumberQ` and `Im == 0` on both
   coefficients so complex coefficients cannot collapse `Sqrt[a^2+b^2]` to zero).
5. `ReplaceRepeated` with `trig_factor_from_sincos` to restore `Tan`/`Sec`/...

**Data structures.** Static rule lists from `trigsimp_init`; tail patterns
(`r___`) on `Plus`/`Times` let each identity fire inside a larger sum/product,
with `ATTR_ORDERLESS` on `Plus`/`Times` driving the permutation search.
Memoized through the active `FactorMemo` via the `builtin_trigfactor` wrapper.

- `Listable`, `Protected`.
- Operates on both circular (`Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc`) and
  hyperbolic (`Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch`) functions.
- Pipeline:
  1. Rewrite reciprocal heads (`Tan`, `Cot`, `Sec`, `Csc`, and their
     hyperbolic analogs) as `Sin`/`Cos`/`Sinh`/`Cosh` ratios so that `Factor`
     sees the full polynomial structure.
  2. Combine into a single rational via `Together`.
  3. Run `Factor` on the resulting rational; trigonometric atoms are treated
     as independent polynomial variables. The `Factor` pass is skipped when
     the post-`Together` form contains more than two distinct squared
     trigonometric atoms (e.g. `Sin[x]^2`, `Cos[x]^2`, `Sinh[y]^2`,
     `Cosh[y]^2` together): on such dense multivariate polynomials Factor's
     trial-division loop stalls without producing a useful factorization, and
     the identity rules in step 4 still match Pythagorean structure that
     survives in the post-`Together` factored form (e.g.
     `(Sin[x]^2 + Cos[x]^2)(Cosh[y]^2 - Sinh[y]^2)` collapses directly to
     `1` via the `Times`-context Pythagorean rules).
  4. Apply identity collapse rules via `ReplaceRepeated`: Pythagorean
     identities (`Sin^2 + Cos^2 -> 1`, `Cosh^2 - Sinh^2 -> 1`, with and
     without arbitrary coefficients), reverse angle-addition
     (`Sin[a]Cos[b] ± Cos[a]Sin[b] -> Sin[a ± b]`,
     `Cos[a]Cos[b] ± Sin[a]Sin[b] -> Cos[a ∓ b]`, and hyperbolic analogs),
     reverse double-angle (`2 Sin Cos -> Sin[2x]`,
     `Cos^2 - Sin^2 -> Cos[2x]`, `Cosh^2 + Sinh^2 -> Cosh[2x]`), and
     factored-form variants such as `(Cos - Sin)(Cos + Sin) -> Cos[2x]`,
     `(Cosh - 1)(Cosh + 1) -> Sinh^2`, and
     `(Cosh - Sinh)(Cosh + Sinh) -> 1` that arise naturally from `Factor`.
  5. Restore `Tan`/`Cot`/`Sec`/`Csc` (and hyperbolic analogs) from the
     `Sin`/`Cos` ratio form so reciprocal heads survive the round-trip.
- Two paths are tried: the primary pipeline (preserves angle-sum structure)
  and a fallback that `TrigExpand`s the argument first (catches
  cancellations that only become visible after the angle-sum is expanded,
  e.g. `Cos[x + y] + Sin[x] Sin[y] -> Cos[x] Cos[y]`). The fallback runs
  only when the primary pipeline leaves the expression unchanged, so
  structurally productive inputs (e.g. `Sin[x + y]^2 + Tan[x + y]`) avoid
  the expensive expanded-rational path. The final result is the smaller of
  the two by leaf count; ties favour the primary pipeline.
- Automatically threads over lists (via `Listable`), as well as equations,
  inequalities (`Equal`, `Unequal`, `Less`, `LessEqual`, `Greater`,
  `GreaterEqual`, `SameQ`, `UnsameQ`), and logic functions (`And`, `Or`,
  `Not`, `Xor`, `Implies`).

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [TrigExpand](../../elementary-functions/TrigExpand/), [Sin](../../elementary-functions/Sin/), [Cos](../../elementary-functions/Cos/), [Tan](../../elementary-functions/Tan/), [Cot](../../elementary-functions/Cot/), [Sec](../../elementary-functions/Sec/), [Csc](../../elementary-functions/Csc/), [Sinh](../../elementary-functions/Sinh/)

- Source: [`src/simp/trigsimp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/trigsimp.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_trigfactor.c`](https://github.com/stblake/mathilda/blob/main/tests/test_trigfactor.c)
