# Element

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Element[x, dom]
    returns True if x is provably an element of the domain dom under the current $Assumptions, False if it is provably not, and stays unevaluated otherwise.
Supported domains: Integers, Rationals, Reals, Algebraics, Complexes, Booleans, Primes, Composites.
Numeric and structural literals decide directly: Element[5, Integers] -> True, Element[5/2, Integers] -> False, Element[1+I, Reals] -> False, Element[2.5, Integers] -> False.
Element consults $Assumptions for symbolic queries, so under Assuming[Element[x, Integers], ...] a query Element[x, Reals] returns True via the Integer => Real lattice.
Element[{x1, ..., xN}, dom] and Element[x1 | ... | xN, dom] are shorthand for the conjunction Element[x1, dom] && ... && Element[xN, dom]: True/False if every component decides, otherwise unevaluated and treated as a joint per-variable fact by Simplify.
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation notes

**Algorithm.** `builtin_element` decides domain membership `Element[x, dom]`,
returning `True`/`False`, or staying unevaluated (NULL) when undetermined.
`Element[{x1,...}, dom]` and `Element[x1|x2|..., dom]` are treated as the
conjunction over components: collapse to `True`/`False` only when every component
decides, otherwise leave the original in place so the assumption framework
(`ctx_walk` in `simp_assume.c`) can still split it into per-variable facts.

For a symbol domain, it reads the current `$Assumptions` (`read_dollar_assumptions`),
builds an `AssumeCtx`, and calls `element_decide`, which returns 1/0/-1. The
decision first checks the fact set directly (`fact_in_domain`), then applies
literal rules per domain: `Integers` (Integer/BigInt yes, integer-valued Real
yes, Rational/Complex no), `Rationals`/`Algebraics`/`Reals`/`Complexes` (with
canonical `Complex` always carrying a non-zero imaginary part, so `Element[I,
Reals]` is False), `Booleans` (only literal `True`/`False`), and
`Primes`/`Composites` (delegating to `PrimeQ`). Provable facts from assumptions
are consulted via `prov_int`/`prov_re`. Unknown or symbolic cases return -1 →
the call stays unevaluated.

**Data structures.** `AssumeCtx` (flat fact array) built per call from
`$Assumptions`; `element_decide` is a `strcmp`-dispatched cascade with no
persistent state.

- `Protected`.
- Supported domains: `Integers`, `Rationals`, `Reals`, `Algebraics`, `Complexes`,

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/simp/simp_builtins.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp_builtins.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
