# DigitSum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DigitSum[n] gives the sum of the decimal digits of the integer n.
DigitSum[n, b] gives the sum of the base-b digits of n.
The sign of n is discarded; DigitSum[0] is 0.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DigitSum[1234]
Out[1]= 10

In[2]:= DigitSum[255, 16]
Out[2]= 30

In[3]:= DigitSum[{1234, 0, 99}]
Out[3]= {10, 0, 18}
```

## Implementation notes

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
