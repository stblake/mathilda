# N

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`N[expr]`**

Gives a machine-precision numerical approximation of expr.

**`N[expr, n]`**

Gives a numerical approximation to n decimal digits. Requires a USE\_MPFR build; without it, a warning is emitted and machine precision is used.

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= N[Sin[3141592653589793238]]
Out[1]= -0.446315

In[2]:= N[Sin[3141592653589793238], 30]
Out[2]= -0.4463151633593201122016036193238
```

### Scope (2)

```mathematica
In[3]:= N[Pi, 100] // N
Out[3]= 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170681

In[4]:= Precision[%]
Out[4]= Infinity
```

### Worked examples (2)

```mathematica
In[5]:= N[1/10^30]
Out[5]= 1e-30

In[6]:= N[10^400/3]
Out[6]= 3.333333333333333e+399
```

### Applications (6)

```mathematica
In[7]:= N[Sqrt[2]]
Out[7]= 1.41421

In[8]:= N[2/7, 15]
Out[8]= 0.2857142857142856

In[9]:= N[Pi, 40]
Out[9]= 3.1415926535897932384626433832795028841971

In[10]:= N[Zeta[3], 40]
Out[10]= 1.2020569031595942853997381615114499907651

In[11]:= N[Gamma[1/3], 35]
Out[11]= 2.67893853470774763365569294097467766

In[12]:= N[EulerGamma, 30]
Out[12]= 0.5772156649015328606065120900823
```

## Algorithm

Mathilda — numeric evaluation implementation.

See numeric.h for the module-level overview and extensibility notes.

This file implements `N[expr]` / `N[expr, prec]`. Phase 1 targets machine-precision IEEE doubles; Phase 2 (gated behind USE_MPFR) adds MPFR arbitrary precision. Phase-2 extension points are marked with an inline "Phase 2" marker so the eventual additions are obvious.

## Implementation notes

**Algorithm.** `builtin_n` parses the optional precision argument (`N[expr]` → machine spec; `N[expr, p]` → `parse_precision_arg`, converting requested decimal digits to MPFR bits via `numeric_digits_to_bits` = `ceil(digits·log2(10))`) into a `NumericSpec`, then calls the recursive `numericalize(expr, spec)`. `numericalize` walks the tree: `EXPR_INTEGER`/`EXPR_BIGINT` become an `EXPR_REAL` (machine mode) or an `EXPR_MPFR` filled at `spec.bits` (MPFR mode); `EXPR_REAL` is promoted to MPFR with zero-padding beyond its 53 exact bits when a higher precision is requested; `EXPR_MPFR` is re-rounded up or down to the target precision, with a guard so a finite MPFR value beyond `DBL_MAX` is kept as a machine-precision MPFR rather than overflowing to ∞. Named constants (`Pi`, `E`, `EulerGamma`, `Catalan`, `GoldenRatio`, `Degree`) are resolved from a registry: a `double` for machine mode, or a dedicated `mpfr_fill` (e.g. `mpfr_const_pi`, with guard digits for derived constants like `Degree = π/180`). Functions are rebuilt with numericalized arguments and re-evaluated, so the actual arithmetic is performed by the MPFR-aware `Plus`/`Times`/`Power`/trig/log kernels.

**Precision propagation.** Precision flows bottom-up through evaluation, not through `N`: each numeric binary op computes its working precision from its operands (`numeric_combined_bits`/`expr_max_mpfr_prec`, with a 53-bit floor) and produces an `EXPR_MPFR` at that precision; `Precision[]`/`Accuracy[]` later report `mpfr_get_prec / log2(10)`. `N` only seeds the leaves at the requested `spec.bits`; mixed-precision results take the minimum-precision contagion from inexact parts. (Without `USE_MPFR`, everything collapses to machine `double`.)

**Data structures.** `NumericSpec { mode, bits }`; arbitrary-precision values are `EXPR_MPFR` wrapping an `mpfr_t`. `N` is registered `LISTABLE | PROTECTED`, so threading over lists happens in the evaluator before the builtin runs.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/), [GoldenRatio](../../mathematical-constants/GoldenRatio/)

- Source: [`src/numeric.c`](https://github.com/stblake/mathilda/blob/main/src/numeric.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_airyai.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airyai.c)
- Tests: [`tests/test_airybi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_airybi.c)
- Tests: [`tests/test_bernoullib.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bernoullib.c)

## Notes & additional examples

### Notes

`N[expr]` gives a machine-precision floating-point value, displayed to about six
significant digits. `N[expr, d]` requests approximately `d` digits of precision,
computed via arbitrary-precision arithmetic (so `N[Pi, 20]` returns the constant
to 20 digits). Exact inputs such as `Sqrt[2]`, `Pi`, `E`, and rationals are
converted to their numeric approximations. Note that machine-precision results
print at the default short width even when more digits are internally available.
