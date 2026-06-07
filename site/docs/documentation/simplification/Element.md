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

- `Protected`.
- Supported domains: `Integers`, `Rationals`, `Reals`, `Algebraics`, `Complexes`,

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
