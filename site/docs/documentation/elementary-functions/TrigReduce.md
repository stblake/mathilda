# TrigReduce

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
TrigReduce[expr]
    rewrites products and powers of trigonometric functions in expr in
    terms of trigonometric functions with combined arguments.
TrigReduce operates on both circular and hyperbolic functions; given a
trigonometric polynomial it typically yields a linear expression
involving trigonometric functions with more complicated arguments
(broadly the inverse of TrigExpand).
TrigReduce automatically threads over lists, equations, inequalities,
and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= TrigReduce[2 Cos[x]^2]
Out[1]= 1 + Cos[2 x]

In[2]:= TrigReduce[2 Sin[x] Cos[y]]
Out[2]= Sin[x + y] + Sin[x - y]

In[3]:= TrigReduce[2 Cosh[x] Cosh[y]]
Out[3]= Cosh[x + y] + Cosh[x - y]

In[4]:= TrigReduce[Sin[a] (Cos[b] - Sin[b]) + Cos[a] (Sin[b] + Cos[b])]
Out[4]= Cos[a + b] + Sin[a + b]

In[5]:= TrigReduce[Tan[x] + Tan[y]]
Out[5]= Sec[x] Sec[y] Sin[x + y]

In[6]:= TrigReduce[Coth[x] + Coth[y]]
Out[6]= Csch[x] Csch[y] Sinh[x + y]

In[7]:= TrigReduce[Sin[x]^4]
Out[7]= 1/8 (3 + Cos[4 x] - 4 Cos[2 x])

In[8]:= TrigReduce[2 Sin[x + y] Cos[x - y]]
Out[8]= Sin[2 x] + Sin[2 y]
```

## Implementation notes

**Algorithm.** `builtin_trigreduce_impl` is the product/power-direction inverse of
TrigExpand: it rewrites products and integer powers of single-argument circular
and hyperbolic trig calls into single trig calls of compound (sum / multiple)
arguments, using the classical product-to-sum and power-reduction identities

```
Sin[a] Cos[b] = (Sin[a+b] + Sin[a-b]) / 2     Cos[a] Cos[b] = (Cos[a+b] + Cos[a-b]) / 2
Sin[a] Sin[b] = (Cos[a-b] - Cos[a+b]) / 2     Sin[x]^2 = (1 - Cos[2x]) / 2
Cos[x]^2 = (1 + Cos[2x]) / 2                   (and the Sinh/Cosh hyperbolic analogues)
```

The pipeline (trig canonicalizer suppressed throughout so the `Sin/Cos`
intermediate forms are not re-collapsed before the rules fire):

1. **To Sin/Cos.** `ReplaceRepeated` with `trig_factor_to_sincos` rewrites
   reciprocal heads (`Tan`/`Cot`/`Sec`/`Csc` and hyperbolic) as `Sin/Cos` ratios
   so the product-to-sum rules can see them.
2. **Iterate to a fixed point** (bounded at 16 iterations): alternate
   `ReplaceRepeated` with `trig_reduce_rules` (the power-reduction and
   product-to-sum identities, with each constructed compound argument wrapped in
   `Expand[...]` so e.g. `Sin[(x+y)-(x-y)]` canonicalizes to `Sin[2y]` before the
   surrounding trig head sees it) and an `Expand` step. The iteration is required
   because `Expand` re-exposes `Cos[2x]^2` terms hidden inside `(1 - Cos[2x])^2/4`
   after a power-reduction pass on `Sin[x]^4`, which the rule then reduces again.
   The 16-iteration cap covers exponents through `Sin[x]^65536`.
3. **Together** to combine over a common denominator so numerators appear as a
   single `Plus` for the collapse rules.
4. **Angle-addition collapse.** `ReplaceRepeated` with `trig_reduce_collapse`:
   coefficient-aware reverse angle-addition (`c. Sin[a]Cos[b] + c. Cos[a]Sin[b]
   :> c Sin[a+b]`, etc.) plus negative-argument cancellation rules guarded by
   `SameQ[Expand[a+b], 0]` that fold `Sin[a-b] + Sin[b-a]` (which the auto-evaluator
   leaves un-reduced) to zero.
5. **From Sin/Cos.** `ReplaceRepeated` with `trig_factor_from_sincos` restores
   `Tan`/`Sec`/`Csc` (and hyperbolic) where the ratio/reciprocal shape survives.
6. **Final canonicalisation:** `Expand` then `Together`, distributing outer scalars
   (`1/2 (2 Cos[a+b] + 2 Sin[a+b])` flattens) while keeping irreducible fractions
   like `(3 - 4 Cos[2x] + Cos[4x])/2` as a single rational.

**Data structures.** Four static rule lists (`trig_factor_to_sincos`,
`trig_reduce_rules`, `trig_reduce_collapse`, `trig_factor_from_sincos`) parsed in
`trigsimp_init`. `Times` is `ATTR_ORDERLESS`, so the matcher commutes factors and
only one direction of each product-to-sum pair needs to be written. Threads over
`List` (via `ATTR_LISTABLE`) and over equation/inequality/logic heads; memoized
through the active `FactorMemo` via the `builtin_trigreduce` wrapper.

- Applies the classical product-to-sum identities (Sin·Cos, Sin·Sin,

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/simp/trigsimp.c`](https://github.com/stblake/mathilda/blob/main/src/simp/trigsimp.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= TrigReduce[Sin[x] Cos[x]]
Out[1]= 1/2 Sin[2 x]
```

A square is linearised through the power-reduction formula:

```mathematica
In[1]:= TrigReduce[Sin[x]^2]
Out[1]= 1/2 (1 - Cos[2 x])
```

Higher powers spread across several harmonics:

```mathematica
In[1]:= TrigReduce[Cos[x]^3]
Out[1]= 1/4 (3 Cos[x] + Cos[3 x])
```

A product of squares reduces to a single fourth harmonic — exactly the
integrand identity behind `Integrate[Sin[x]^2 Cos[x]^2, x]`:

```mathematica
In[1]:= TrigReduce[Sin[x]^2 Cos[x]^2]
Out[1]= 1/8 (1 - Cos[4 x])
```

The product-to-sum identity for two sines appears directly:

```mathematica
In[1]:= TrigReduce[2 Sin[x] Sin[y]]
Out[1]= -Cos[x + y] + Cos[x - y]
```
