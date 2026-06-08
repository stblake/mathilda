# Round

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Round[x]
    rounds x to the nearest integer, breaking ties to the nearest even
    integer (banker's rounding).
Round[x, a]
    rounds x to the nearest multiple of a.
Round is Listable. Exact inputs return exact integers; Real / MPFR
inputs round at the input precision.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_round` is `do_piecewise(res, OP_ROUND, "Round", allow_2_args=true)`, sharing the rounding kernel with `Floor`/`Ceiling`/`IntegerPart`/`FractionalPart`. Round uses **banker's rounding** (round-half-to-even) everywhere. For machine `double`/`EXPR_REAL` values it calls `round_half_even`, which floors and inspects the fractional residue, breaking exact `.5` ties toward the even integer. For exact rationals it works entirely in GMP: it computes `floor((2·num + den) / (2·den))` (round-half-up) and corrects the exact-half tie (`rem == 0` and `q` odd) by decrementing `q` toward the even neighbour. MPFR values round via `mpfr_rint(..., MPFR_RNDN)`. Complex values round real and imaginary parts independently. An exact non-rational numeric quantity (e.g. `10000000 Pi`) is resolved through `do_piecewise_numeric_exact`, which numericalises at growing precision until the integer result is unambiguous.

The two-argument form `Round[x, a]` is rewritten as `a · Round[x/a]` (rounding to the nearest multiple of `a`). Symbolic arguments return `NULL`. Attributes include `ATTR_LISTABLE | ATTR_NUMERICFUNCTION`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/piecewise.c`](https://github.com/stblake/mathilda/blob/main/src/piecewise.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Round[7/2]
Out[1]= 4

In[2]:= Round[5/2]
Out[2]= 2

In[3]:= Round[17, 5]
Out[3]= 15

In[4]:= Round[{1.4, 2.5, 3.6}]
Out[4]= {1, 2, 4}
```

### Notes

`Round` breaks ties to the nearest even integer (banker's rounding), so `Round[5/2]` and `Round[2.5]` both give 2. `Round[x, a]` rounds to the nearest multiple of `a`; Round is Listable.
