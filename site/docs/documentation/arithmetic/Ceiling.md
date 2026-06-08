# Ceiling

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Ceiling[x]
    gives the smallest integer greater than or equal to x.
Ceiling[x, a]
    gives the smallest multiple of a greater than or equal to x.
Ceiling is Listable. Exact inputs return exact integers; Real / MPFR
inputs are rounded toward +Infinity at the input precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_ceiling` (`src/piecewise.c`) is a thin wrapper over `do_piecewise(res, OP_CEILING, ...)`, shared with Floor/Round/IntegerPart/FractionalPart. The 1-arg form dispatches to `do_piecewise_1`, which handles each exact kind directly: integers/bigints pass through; `EXPR_REAL` uses C `ceil`; `EXPR_MPFR` uses `mpfr_ceil`; exact `Rational[p,q]` (including bigint components) uses GMP `mpz_cdiv_q` for an exact integer result; `Complex[re,im]` recurses componentwise. For an exact-but-symbolic real numeric that the leaf branches cannot resolve, `do_piecewise_numeric_exact` numericalizes to MPFR at doubling precision (256 up to 2^16 bits) and accepts the ceiling only once two successive precisions agree — an interval-style certification that avoids mis-rounding values near an integer boundary, returning `NULL` rather than a wrong answer if it cannot converge. The 2-arg `Ceiling[x, a]` rewrites to `a * Ceiling[x/a]`.

**Data structures.** `Expr*` leaves plus GMP `mpz_t`/`mpq` arithmetic and MPFR `mpfr_t` for the precision-escalation certification.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/piecewise.c`](https://github.com/stblake/mathilda/blob/main/src/piecewise.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Ceiling[7/2]
Out[1]= 4

In[2]:= Ceiling[-2.7]
Out[2]= -2

In[3]:= Ceiling[17, 5]
Out[3]= 20
```

### Notes

`Ceiling[x]` rounds toward `+Infinity`; the two-argument `Ceiling[x, a]` gives the smallest multiple of `a` at least `x`. Exact inputs return exact integers.
