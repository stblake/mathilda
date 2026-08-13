# TrigExpand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TrigExpand[expr]`**

expands out trigonometric functions in expr. TrigExpand operates on both circular and hyperbolic functions. TrigExpand splits up sums and integer multiples that appear in arguments of trigonometric functions, and then expands out products of trigonometric functions into sums of powers, using trigonometric identities when possible. TrigExpand automatically threads over lists, as well as equations, inequalities, and logic functions.

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= TrigExpand[Sin[2 x]]
Out[1]= 2 Cos[x] Sin[x]

In[2]:= TrigExpand[Sin[x + y]]
Out[2]= Cos[x] Sin[y] + Sin[x] Cos[y]

In[3]:= TrigExpand[Sin[3 x]]
Out[3]= -Sin[x]^3 + 3 Cos[x]^2 Sin[x]

In[4]:= TrigExpand[Cos[x + y + z]]
Out[4]= Cos[x] Cos[y] Cos[z] - Cos[x] Sin[y] Sin[z] - Sin[x] Cos[y] Sin[z] - Sin[x] Sin[y] Cos[z]

In[5]:= TrigExpand[Sinh[4 x]]
Out[5]= 4 Cosh[x] Sinh[x]^3 + 4 Cosh[x]^3 Sinh[x]

In[6]:= TrigExpand[Cosh[x - y]]
Out[6]= Cosh[x] Cosh[y] - Sinh[x] Sinh[y]

In[7]:= TrigExpand[Tanh[2 t]]
Out[7]= 2 Cosh[t] Sinh[t] Sech[2 t]
```

### Applications (4)

```mathematica
In[8]:= TrigExpand[Sin[a + b]]
Out[8]= Cos[a] Sin[b] + Sin[a] Cos[b]

In[9]:= TrigExpand[Sin[3 x]]
Out[9]= -Sin[x]^3 + 3 Cos[x]^2 Sin[x]

In[10]:= TrigExpand[Cos[2 x]]
Out[10]= Cos[x]^2 - Sin[x]^2

In[11]:= TrigExpand[Sinh[x + y]]
Out[11]= Cosh[x] Sinh[y] + Sinh[x] Cosh[y]
```

## Implementation notes

**Algorithm.** `builtin_trigexpand_impl` expands trig/hyperbolic functions of
sums and integer multiples into products and powers of single-argument trig
calls. The pipeline (with the trig canonicalizer suppressed throughout):

1. **Angle-addition + multiple-angle.** `ReplaceRepeated` with
   `trig_expand_rules`. The binary angle-addition forms
   `Sin[x_ + y__] :> Sin[x] Cos[Plus[y]] + Cos[x] Sin[Plus[y]]` (and `Cos`,
   `Sinh`, `Cosh`) recurse over the rest of the summands; the multiple-angle
   forms `Sin[n_Integer x_] /; n>1 :> Sin[(n-1)x] Cos[x] + Cos[(n-1)x] Sin[x]`
   reduce integer multiples down to `Sin[x]`/`Cos[x]`. Reciprocal heads
   (`Tan`/`Cot`/`Sec`/`Csc` and hyperbolic analogues) are rewritten as `Sin/Cos`
   ratios so the base rules apply. A large block of inverse-trig composition
   rules (`Cos[ArcSin[x]] :> Sqrt[1-x^2]`, etc.) is included.
2. **Expand** to distribute products of sums into a flat monomial sum.
3. **Pythagorean collapse.** If the expanded form has a denominator
   (`has_reciprocal_power`), no Pythagorean-eligible squared pair
   (`input_has_pythag_pair`), or more than `TRIG_FACTOR_ATOM_THRESHOLD` distinct
   squared trig atoms, only the direct-sum rules `trig_expand_pythag` are applied
   (`ReplaceRepeated`). Otherwise it first runs polynomial `Factor` — which turns
   `Sin[nx]^2 + Cos[nx]^2` into `(Sin[x]^2+Cos[x]^2)^n` — then applies
   `trig_expand_pythag` to collapse `(Sin^2+Cos^2)^n -> 1` (and the negated-sign
   and hyperbolic variants).
4. **Re-Expand** to restore the canonical monomial form.

**Data structures.** `trig_expand_rules` and `trig_expand_pythag` are static
`parse_expression`'d rule lists built in `trigsimp_init`. Threads over `List`
(via `ATTR_LISTABLE`) and over equations/inequalities/logic heads
(`trigexpand_threads_over`). Memoized through the active `FactorMemo` by the
`builtin_trigexpand` wrapper (`trig_memo_call`).

- `Listable`, `Protected`.
- Operates on both circular (`Sin`, `Cos`, `Tan`, `Cot`, `Sec`, `Csc`) and
  hyperbolic (`Sinh`, `Cosh`, `Tanh`, `Coth`, `Sech`, `Csch`) functions.
- Applies angle-addition formulas to `Sin[a + b + …]`, `Cos[a + b + …]`,
  `Sinh[a + b + …]`, `Cosh[a + b + …]` to a fixed point.
- Applies multiple-angle reductions to `Sin[n x]`, `Cos[n x]`, `Sinh[n x]`,
  `Cosh[n x]` for integer `n`, recursively reducing to `Sin[x]` / `Cos[x]` /
  `Sinh[x]` / `Cosh[x]`.
- `Tan`, `Cot`, `Sec`, `Csc` (and `Tanh`, `Coth`, `Sech`, `Csch`) with sum or
  integer-multiple arguments are rewritten as ratios of `Sin`/`Cos`
  (resp. `Sinh`/`Cosh`) and then expanded.
- Distributes products over sums via `Expand` so the result is a flat sum of
  monomials.
- Applies the Pythagorean identities `Sin[x]^2 + Cos[x]^2 -> 1` and
  `Cosh[x]^2 - Sinh[x]^2 -> 1` as a final reduction pass, including powers of
  both identities for any integer `n >= 1`:
    - `Sin[n x]^2 + Cos[n x]^2` expands to `(Sin[x]^2 + Cos[x]^2)^n` and
      collapses to `1` via a Factor-based reduction.
    - `Cosh[n x]^2 - Sinh[n x]^2` factors as
      `(Cosh[x] + Sinh[x])^n (Cosh[x] - Sinh[x])^n` and collapses to `1`.
  Negated and scalar-weighted forms (e.g. `-Sin[n x]^2 - Cos[n x]^2`,
  `-5 (Sin[n x]^2 + Cos[n x]^2)`, `Sinh[n x]^2 - Cosh[n x]^2`) collapse to the
  expected signed constant — the Pythagorean rules match both possible signs
  that `Factor` may emerge with and allow an arbitrary remainder of factors in
  the surrounding `Times`. Expressions that contain a denominator (any
  `Power[_, negative_Integer]` subterm) skip the Factor pass so that canonical
  forms such as `(2 Cos[x] Sin[x])/(Cos[x]^2 - Sin[x]^2)` are preserved.
  Inputs without a Pythagorean-eligible squared structure (no pair
  `Sin[a]^k`/`Cos[a]^k` or `Sinh[a]^k`/`Cosh[a]^k` with the same argument and
  `k >= 2`) likewise skip the Factor pass; the multivariate polynomials that
  multi-angle expansions such as `TrigExpand[Sin[2 x + 3 y]]` produce would
  otherwise make `Factor` prohibitively slow without yielding any collapse.
  The Factor pass is also skipped when the expanded form contains more than
  two distinct squared trigonometric atoms (e.g. `Cos[x]^2`, `Sin[x]^2`,
  `Cos[y]^2`, `Sin[y]^2` together): even if a Pythagorean pair is structurally
  present, `Factor` on the resulting dense multivariate polynomial stalls
  without producing a useful collapse.
- Automatically threads over lists (via `Listable`), as well as equations,
  inequalities (`Equal`, `Unequal`, `Less`, `LessEqual`, `Greater`,
  `GreaterEqual`, `SameQ`, `UnsameQ`), and logic functions (`And`, `Or`,
  `Not`, `Xor`, `Implies`).

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Sin](../../elementary-functions/Sin/), [Cos](../../elementary-functions/Cos/), [Tan](../../elementary-functions/Tan/), [Cot](../../elementary-functions/Cot/), [Sec](../../elementary-functions/Sec/), [Csc](../../elementary-functions/Csc/), [Sinh](../../elementary-functions/Sinh/), [Cosh](../../elementary-functions/Cosh/)

- Source: [`src/simp/trigsimp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/trigsimp.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
- Tests: [`tests/test_intrischnorman.c`](https://github.com/stblake/mathilda/blob/main/tests/test_intrischnorman.c)
- Tests: [`tests/test_trigexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_trigexpand.c)
- Tests: [`tests/test_trigfactor.c`](https://github.com/stblake/mathilda/blob/main/tests/test_trigfactor.c)
