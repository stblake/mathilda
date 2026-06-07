# FromContinuedFraction

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FromContinuedFraction[{a1, a2, ..., an}]
    reconstructs a1 + 1/(a2 + 1/(a3 + ... + 1/an)). The terms may be
    symbolic; the result is the convergent in nested (un-expanded) form.
FromContinuedFraction[{a1, ..., am, {b1, ..., bk}}]
    gives the exact quadratic irrational whose continued-fraction terms
    begin with the ai then cycle through the bi forever; all ai and bi
    must be integers. FromContinuedFraction[{}] is 0. It is the inverse
    of ContinuedFraction.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FromContinuedFraction[{2, 1, 3, 4}]
Out[1]= 47/17

In[2]:= FromContinuedFraction[{a, b, c, d}]
Out[2]= (1 + a b + (a + (1 + a b) c) d)/(b + (1 + b c) d)

In[3]:= FromContinuedFraction[{8, {2, 2, 1, 7, 1, 2, 2, 16}}]
Out[3]= Sqrt[71]

In[4]:= FromContinuedFraction[{{1, 2, 3, 4}}]
Out[4]= 1/15 (9 + 2 Sqrt[39])

In[5]:= FromContinuedFraction[ContinuedFraction[Pi, 3]]
Out[5]= 333/106
```

## Implementation notes

- `Protected` (not `Listable` — the argument is the whole term list).
- The `ai` of the finite form may be **symbolic**; the result is the convergent

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
