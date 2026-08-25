# FractionalPart

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FractionalPart[x]`**

gives the fractional part x - IntegerPart\[x\], carrying the sign of x, so that FractionalPart\[2.7\] is 0.7 and FractionalPart\[-2.7\] is -0.7.

<details>
<summary>Notes</summary>

FractionalPart is Listable and preserves the precision of x. Exact inputs stay exact; symbolic inputs stay unevaluated.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= FractionalPart[2.7]
Out[1]= 0.7

In[2]:= FractionalPart[-2.7]
Out[2]= -0.7

In[3]:= FractionalPart[7/2]
Out[3]= 1/2

In[4]:= IntegerPart[-2.7] + FractionalPart[-2.7]
Out[4]= -2.7
```

## Implementation notes

`builtin_fractionalpart` computes `x - IntegerPart[x]` through the shared `do_piecewise(res, OP_FRACPART, ...)` kernel, keeping the sign of `x` and the precision of the input: `EXPR_REAL` returns `v - trunc(v)`, `EXPR_MPFR` subtracts `mpfr_trunc` at full precision, and exact rationals return an exact `Rational`. Registered `PROTECTED | NUMERICFUNCTION | LISTABLE`; a quantity with no monotone reduction (e.g. `FractionalPart[10^7 3^(2/3)]`) is left symbolic.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/piecewise.c`](https://github.com/stblake/mathilda/blob/main/src/piecewise.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)

## Notes & additional examples

### Notes

`FractionalPart[x]` is `x - IntegerPart[x]`, so it carries the **sign of `x`**:
`FractionalPart[-2.7] = -0.7`, not `0.3`. It preserves the input's precision and
keeps exact inputs exact (`FractionalPart[7/2] = 1/2`). `FractionalPart` is
`Listable`, and reconstructs the number with `IntegerPart`.
