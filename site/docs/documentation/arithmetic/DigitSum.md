# DigitSum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DigitSum[n] gives the sum of the decimal digits of the integer n.`**

**`DigitSum[n, b] gives the sum of the base-b digits of n.`**

<details>
<summary>Notes</summary>

The sign of n is discarded; DigitSum\[0\] is 0.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

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

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_digit_sum.c`](https://github.com/stblake/mathilda/blob/main/tests/test_digit_sum.c)
