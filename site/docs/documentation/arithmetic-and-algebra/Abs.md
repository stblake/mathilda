# Abs

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Abs[z] gives the absolute value (modulus) of numeric z, Sqrt[Re[z]^2 + Im[z]^2] for complex z.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_abs` handles each numeric kind directly: `mpz_abs` for `EXPR_BIGINT`, sign flips for `EXPR_INTEGER`/`EXPR_REAL`/`EXPR_MPFR` (via `mpfr_abs`), and `|n|/d` for rationals. For a `Complex[re, im]` literal (or an expression `complex_decompose` splits into numeric real/imag parts) it builds the symbolic modulus `Power[Plus[re^2, im^2], 1/2]`; when MPFR components are present it folds directly through `mpfr_hypot` at the combined working precision instead, which is also numerically stable across disparate magnitudes. Symbolic arguments return `NULL`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
