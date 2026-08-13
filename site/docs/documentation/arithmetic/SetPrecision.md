# SetPrecision

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SetPrecision[x, n]`**

Returns an expression equivalent to x with numeric values re-rounded or promoted to n decimal digits of precision. Requires a USE\_MPFR build for n \> MachinePrecision.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= SetPrecision[1.5, 30]
Out[1]= 1.5

In[2]:= Precision[SetPrecision[1.5, 30]]
Out[2]= 30.103

In[3]:= SetPrecision[Pi, 50]
Out[3]= 3.1415926535897932384626433832795028841971693993751

In[4]:= SetPrecision[1/3, 40]
Out[4]= 0.33333333333333333333333333333333333333332
```

## Implementation notes

`builtin_set_precision` is a two-argument wrapper: it parses the precision argument into a `NumericSpec` via `parse_prec_arg` (accepting an integer/real digit count or `MachinePrecision`) and drives `numericalize(value, spec)` — the same engine `N` uses — to re-represent the value to the requested number of significant digits (an `EXPR_MPFR` at the corresponding bit width when MPFR is built, otherwise machine `double`). Returns `NULL` if the precision argument is not a valid spec.

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/precision.c`](https://github.com/stblake/mathilda/blob/main/src/precision.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_mateigen_direct.c`](https://github.com/stblake/mathilda/blob/main/tests/test_mateigen_direct.c)
- Tests: [`tests/test_matinv_methods.c`](https://github.com/stblake/mathilda/blob/main/tests/test_matinv_methods.c)
- Tests: [`tests/test_matsol_methods.c`](https://github.com/stblake/mathilda/blob/main/tests/test_matsol_methods.c)
- Tests: [`tests/test_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric.c)

## Notes & additional examples

### Notes

`SetPrecision[x, n]` returns a value equal to `x` but carrying `n` digits of precision; the printed form may look unchanged while the internal precision is raised (confirm with `Precision`). Padding extra digits onto a machine-precision number introduces meaningless trailing bits, so only widen precision when the original value genuinely has them.
