# SetAccuracy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SetAccuracy[x, n]`**

Returns an expression equivalent to x with numeric values re-rounded or promoted to n decimal digits of accuracy. Requires a USE\_MPFR build for high-accuracy outputs.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= SetAccuracy[1.5, 30]
Out[1]= 1.5

In[2]:= Accuracy[SetAccuracy[1.5, 30]]
Out[2]= 30.2279

In[3]:= SetAccuracy[Pi, 30]
Out[3]= 3.141592653589793238462643383279

In[4]:= SetAccuracy[123.456, 20]
Out[4]= 123.45600000000000306954
```

## Implementation notes

`builtin_set_accuracy` re-expresses a value to a target *accuracy* (digits after the decimal point) by converting accuracy to precision. It extracts the numeric accuracy `n` (integer/real/rational, or `MachinePrecision` which short-circuits to a machine-spec `numericalize`), then computes the required significant digits as `digits = n + log10(|x|)` using `expr_log10_abs`, floored at 1. It builds a `NumericSpec` (MPFR bits via `numeric_digits_to_bits(digits)`, or machine spec without MPFR) and calls `numericalize`. This is the standard "accuracy = digits past the point" approximation, not full significance-arithmetic semantics. Non-positive accuracy or unrecognised argument types return `NULL`.

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/precision.c`](https://github.com/stblake/mathilda/blob/main/src/precision.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric.c)

## Notes & additional examples

### Notes

`SetAccuracy[x, n]` returns a value equal to `x` with `n` digits of accuracy (digits past the decimal point); use `Accuracy` to confirm, since the printed form often looks unchanged. It is the absolute-magnitude counterpart to `SetPrecision`.
