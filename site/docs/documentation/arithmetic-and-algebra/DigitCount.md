# DigitCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DigitCount[n] gives a list of the counts of digits 1, 2, ..., 9, 0 in the base-10 representation of n.
DigitCount[n, b] gives a list of the counts of digits 1, 2, ..., b-1, 0 in the base-b representation of n.
DigitCount[n, b, d] gives the number of d digits in the base-b representation of n.
The sign of n is discarded; DigitCount[0] is a list of zeros.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DigitCount[2147, 2, 1]
Out[1]= 5

In[2]:= DigitCount[2147, 2]
Out[2]= {5, 7}

In[3]:= DigitCount[100!]
Out[3]= {15, 19, 10, 10, 14, 19, 7, 14, 20, 30}
```

## Implementation notes

- `Protected` (intentionally not `Listable`). `DigitCount[{1,2,3}]` is left

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
