# FindIntegerNullVector

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FindIntegerNullVector[{x1, ..., xn}]
    finds integers {a1, ..., an}, not all zero, with a1 x1 + ... + an xn == 0 (PSLQ / integer-relation detection).
FindIntegerNullVector[{x1, ..., xn}, d]
    restricts the search to relations of norm <= d.
The xi may be real or complex, exact or inexact; for complex xi the ai are Gaussian integers.  Exact relations are validated with PossibleZeroQ; for inexact xi the relation holds to the precision of the input.  When no relation is found the call is returned unevaluated.
Options:
    WorkingPrecision    Automatic, or a digit count for the search.
    ZeroTest            Automatic, or a function applied to the residual.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FindIntegerNullVector[{Log[2], Log[4]}]
Out[1]= {-2, 1}

In[2]:= FindIntegerNullVector[{Pi, ArcTan[1/5], ArcTan[1/239]}]
Out[2]= {1, -16, 4}

In[3]:= a = Sqrt[2] + 3^(1/3); FindIntegerNullVector[a^Range[0, 6]]
Out[3]= {1, -36, 12, -6, -6, 0, 1}

In[4]:= FindIntegerNullVector[{1, 2 I + Sqrt[3], (2 I + Sqrt[3])^2}]
Out[4]= {-7, -4*I, 1}

In[5]:= FindIntegerNullVector[{E, Pi}, 1000000]
Out[5]= FindIntegerNullVector[{E, Pi}, 1000000]
```

## Implementation notes

- `Protected`.
- The `xi` may be **real or complex**, **exact or inexact**. For complex

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
