# N

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
N[expr]
    Gives a machine-precision numerical approximation of expr.
N[expr, n]
    Gives a numerical approximation to n decimal digits. Requires
    a USE_MPFR build; without it, a warning is emitted and machine
    precision is used.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= N[Pi, 100] // N
Out[1]= 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170681

In[2]:= Precision[%]
Out[2]= 100.243
```

## Implementation notes

**Algorithm.** `builtin_n` parses the optional precision argument (`N[expr]` → machine spec; `N[expr, p]` → `parse_precision_arg`, converting requested decimal digits to MPFR bits via `numeric_digits_to_bits` = `ceil(digits·log2(10))`) into a `NumericSpec`, then calls the recursive `numericalize(expr, spec)`. `numericalize` walks the tree: `EXPR_INTEGER`/`EXPR_BIGINT` become an `EXPR_REAL` (machine mode) or an `EXPR_MPFR` filled at `spec.bits` (MPFR mode); `EXPR_REAL` is promoted to MPFR with zero-padding beyond its 53 exact bits when a higher precision is requested; `EXPR_MPFR` is re-rounded up or down to the target precision, with a guard so a finite MPFR value beyond `DBL_MAX` is kept as a machine-precision MPFR rather than overflowing to ∞. Named constants (`Pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio`, `Degree`) are resolved from a registry: a `double` for machine mode, or a dedicated `mpfr_fill` (e.g. `mpfr_const_pi`, with guard digits for derived constants like `Degree = π/180`). Functions are rebuilt with numericalized arguments and re-evaluated, so the actual arithmetic is performed by the MPFR-aware `Plus`/`Times`/`Power`/trig/log kernels.

**Precision propagation.** Precision flows bottom-up through evaluation, not through `N`: each numeric binary op computes its working precision from its operands (`numeric_combined_bits`/`expr_max_mpfr_prec`, with a 53-bit floor) and produces an `EXPR_MPFR` at that precision; `Precision[]`/`Accuracy[]` later report `mpfr_get_prec / log2(10)`. `N` only seeds the leaves at the requested `spec.bits`; mixed-precision results take the minimum-precision contagion from inexact parts. (Without `USE_MPFR`, everything collapses to machine `double`.)

**Data structures.** `NumericSpec { mode, bits }`; arbitrary-precision values are `EXPR_MPFR` wrapping an `mpfr_t`. `N` is registered `LISTABLE | PROTECTED`, so threading over lists happens in the evaluator before the builtin runs.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/numeric.c`](https://github.com/stblake/mathilda/blob/main/src/numeric.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= N[Pi, 20]
Out[1]= 3.14159265358979323846
```

```mathematica
In[1]:= N[Sqrt[2]]
Out[1]= 1.41421
```

```mathematica
In[1]:= N[2/7, 15]
Out[1]= 0.2857142857142856
```

```mathematica
In[1]:= N[E]
Out[1]= 2.71828
```

### Notes

`N[expr]` gives a machine-precision floating-point value, displayed to about six
significant digits. `N[expr, d]` requests approximately `d` digits of precision,
computed via arbitrary-precision arithmetic (so `N[Pi, 20]` returns the constant
to 20 digits). Exact inputs such as `Sqrt[2]`, `Pi`, `E`, and rationals are
converted to their numeric approximations. Note that machine-precision results
print at the default short width even when more digits are internally available.
