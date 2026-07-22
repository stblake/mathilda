# Precision

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Precision[x]
    Returns the number of decimal digits of precision in x.
    Exact numbers return Infinity; machine-precision reals return
    the symbol MachinePrecision; MPFR values return their declared
    precision in decimal digits.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_precision` (1-arg) delegates to the recursive `precision_of`. Exact quantities — `EXPR_INTEGER`, `EXPR_BIGINT`, exact `Rational`, exact numeric symbols — return `Infinity`. `EXPR_REAL` returns the symbol `MachinePrecision`. `EXPR_MPFR` returns its decimal precision computed from the stored bit precision (`mpfr_get_prec / log2(10)`) as an `EXPR_REAL`. For composite expressions, `Complex[re, im]` and general function heads take the minimum precision across parts via `precision_min` (which treats `MachinePrecision` as the constant `NUMERIC_MACHINE_PRECISION_DIGITS ≈ 15.95` when comparing against explicit MPFR digit counts). This follows the rule that an expression is only as precise as its least-precise inexact part; the precision-tracking unit conversion (`LOG2_10`) is shared with `numeric.c`.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/precision.c`](https://github.com/stblake/mathilda/blob/main/src/precision.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Precision[N[Pi, 30]]
Out[1]= 30.103

In[2]:= Precision[1]
Out[2]= Infinity

In[3]:= Precision[1.5]
Out[3]= MachinePrecision
```

Arithmetic is precision-contagious: a sum is no more precise than its least precise operand, so adding a 30-digit number to a 50-digit number yields about 30 digits:

```mathematica
In[1]:= Precision[N[Pi, 50] + N[E, 30]]
Out[1]= 30.103
```

Squaring a 100-digit square root *gains* a fraction of a digit, reflecting the conditioning of the operation:

```mathematica
In[1]:= Precision[N[Sqrt[2], 100]^2]
Out[1]= 100.243
```

### Notes

`Precision` gives the number of significant decimal digits in a number. Exact quantities (integers, rationals, symbols) have `Infinity` precision, machine-precision reals report `MachinePrecision`, and arbitrary-precision reals report their actual digit count.
