# DifferenceDelta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DifferenceDelta[f, i] gives the forward difference (f /. i -> i+1) - f, the discrete analogue of D. It is the left inverse of indefinite Sum.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DifferenceDelta[i^2, i]
Out[1]= 1 + 2 i

In[2]:= DifferenceDelta[Sum[k k!, k], k]
Out[2]= -Factorial[k] + Factorial[1 + k]
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/sum/sum_gosper.c`](https://github.com/stblake/mathilda/blob/main/src/sum/sum_gosper.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)
