# SetAccuracy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SetAccuracy[x, n]
    Returns an expression equivalent to x with numeric values
    re-rounded or promoted to n decimal digits of accuracy.
    Requires a USE_MPFR build for high-accuracy outputs.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_set_accuracy` re-expresses a value to a target *accuracy* (digits after the decimal point) by converting accuracy to precision. It extracts the numeric accuracy `n` (integer/real/rational, or `MachinePrecision` which short-circuits to a machine-spec `numericalize`), then computes the required significant digits as `digits = n + log10(|x|)` using `expr_log10_abs`, floored at 1. It builds a `NumericSpec` (MPFR bits via `numeric_digits_to_bits(digits)`, or machine spec without MPFR) and calls `numericalize`. This is the standard "accuracy = digits past the point" approximation, not full significance-arithmetic semantics. Non-positive accuracy or unrecognised argument types return `NULL`.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/precision.c`](https://github.com/stblake/mathilda/blob/main/src/precision.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= SetAccuracy[1.5, 30]
Out[1]= 1.5

In[2]:= Accuracy[SetAccuracy[1.5, 30]]
Out[2]= 30.2279
```

Promoting an exact constant to a fixed accuracy yields its high-accuracy decimal
expansion — here 30 digits past the point for `Pi`:

```mathematica
In[1]:= SetAccuracy[Pi, 30]
Out[1]= 3.141592653589793238462643383279
```

Forcing accuracy *beyond* the genuine precision of a machine number exposes the
binary round-off in its tail (the digits past `123.456` are noise):

```mathematica
In[1]:= SetAccuracy[123.456, 20]
Out[1]= 123.45600000000000306954
```

### Notes

`SetAccuracy[x, n]` returns a value equal to `x` with `n` digits of accuracy (digits past the decimal point); use `Accuracy` to confirm, since the printed form often looks unchanged. It is the absolute-magnitude counterpart to `SetPrecision`.
