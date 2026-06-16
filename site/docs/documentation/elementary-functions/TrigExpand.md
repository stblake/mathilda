# TrigExpand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TrigExpand[expr]
    expands out trigonometric functions in expr.
    TrigExpand operates on both circular and hyperbolic functions.
    TrigExpand splits up sums and integer multiples that appear in arguments of
    trigonometric functions, and then expands out products of trigonometric
    functions into sums of powers, using trigonometric identities when possible.
    TrigExpand automatically threads over lists, as well as equations,
    inequalities, and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= TrigExpand[Sin[2 x]]
Out[1]= 2 Cos[x] Sin[x]

In[2]:= TrigExpand[Sin[x + y]]
Out[2]= Cos[x] Sin[y] + Sin[x] Cos[y]

In[3]:= TrigExpand[Sin[3 x]]
Out[3]= -Sin[x]^3 + 3 Cos[x]^2 Sin[x]

In[4]:= TrigExpand[Cos[x + y + z]]
Out[4]= Cos[x] Cos[y] Cos[z] - Cos[x] Sin[y] Sin[z] - Sin[x] Cos[y] Sin[z] - Sin[x] Sin[y] Cos[z]

In[5]:= TrigExpand[Sin[x]^2 + Cos[x]^2]
Out[5]= 1

In[6]:= TrigExpand[Sinh[4 x]]
Out[6]= 4 Cosh[x] Sinh[x]^3 + 4 Cosh[x]^3 Sinh[x]

In[7]:= TrigExpand[Cosh[x - y]]
Out[7]= Cosh[x] Cosh[y] - Sinh[x] Sinh[y]

In[8]:= TrigExpand[Tanh[2 t]]
Out[8]= 2 Cosh[t] Sinh[t] Sech[2 t]
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

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/simp/trigsimp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/trigsimp.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= TrigExpand[Sin[a + b]]
Out[1]= Cos[a] Sin[b] + Sin[a] Cos[b]
```

Integer multiples in the argument are expanded via the multiple-angle formulas:

```mathematica
In[1]:= TrigExpand[Sin[3 x]]
Out[1]= -Sin[x]^3 + 3 Cos[x]^2 Sin[x]
```

The double-angle identity for the cosine drops out automatically:

```mathematica
In[1]:= TrigExpand[Cos[2 x]]
Out[1]= Cos[x]^2 - Sin[x]^2
```

`TrigExpand` works on hyperbolic functions as well, giving the addition formula
for `sinh`:

```mathematica
In[1]:= TrigExpand[Sinh[x + y]]
Out[1]= Cosh[x] Sinh[y] + Sinh[x] Cosh[y]
```
