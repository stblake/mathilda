# ContinuedFraction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ContinuedFraction[x, n]
    gives a list of the first n terms in the continued-fraction
    representation of x.
ContinuedFraction[x]
    gives all terms determinable from the precision of x.
The list {a1, a2, a3, ...} corresponds to a1 + 1/(a2 + 1/(a3 + ...)).
Exact rationals give a finite (canonical, last term >= 2) expansion.
For Sqrt[d] with d a non-square integer the no-count form returns
{a1, ..., {b1, ...}}, the bracketed block repeating cyclically. Inexact
Real / MPFR inputs yield terms only as far as the precision determines
them. ContinuedFraction is Listable.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ContinuedFraction[47/17]
Out[1]= {2, 1, 3, 4}

In[2]:= ContinuedFraction[Sqrt[13]]
Out[2]= {3, {1, 1, 1, 1, 6}}

In[3]:= ContinuedFraction[Pi, 20]
Out[3]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2}

In[4]:= ContinuedFraction[N[Pi]]
Out[4]= {3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14}

In[5]:= ContinuedFraction[Exp[Pi Sqrt[163]], 10]
Out[5]= {262537412640768743, 1, 1333462407511, 1, 8, 1, 1, 5, 1, 4}
```

## Implementation notes

- `Protected`, `Listable`.
- **Exact rationals** (Integer / BigInt / Rational) use the Euclidean

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
