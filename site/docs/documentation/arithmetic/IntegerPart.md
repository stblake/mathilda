# IntegerPart

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IntegerPart[x]`**

gives the integer part of x, truncated toward zero, so that IntegerPart\[2.7\] is 2 and IntegerPart\[-2.7\] is -2.

**`IntegerPart[x] + FractionalPart[x] == x.`**

<details>
<summary>Notes</summary>

IntegerPart is Listable. Exact (Integer / BigInt / Rational) inputs return exact integers; Real / MPFR inputs are truncated at the input precision; symbolic inputs stay unevaluated. It satisfies

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= IntegerPart[2.7]
Out[1]= 2

In[2]:= IntegerPart[-2.7]
Out[2]= -2

In[3]:= IntegerPart[7/2]
Out[3]= 3

In[4]:= IntegerPart[{2.7, -2.7, 7/2, -7/2}]
Out[4]= {2, -2, 3, -3}
```

## Implementation notes

`builtin_integerpart` dispatches through the shared `do_piecewise(res, OP_INTPART, ...)` kernel, truncating each numeric type *toward zero*: Integer/BigInt pass through; rationals use `mpz_tdiv_q` on numerator and denominator; `EXPR_REAL` uses C `trunc()`; `EXPR_MPFR` uses `mpfr_trunc` then `mpfr_get_z` into an `mpz_t`, so arbitrarily large values never silently overflow `int64_t`. Registered `PROTECTED | NUMERICFUNCTION | LISTABLE`. It satisfies `IntegerPart[x] + FractionalPart[x] == x`; non-numeric arguments stay symbolic.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/piecewise.c`](https://github.com/stblake/mathilda/blob/main/src/piecewise.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_interval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_interval.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`IntegerPart[x]` truncates **toward zero**, so `IntegerPart[-2.7] = -2` — unlike
`Floor`, which rounds toward `-Infinity` and would give `-3`. Exact inputs return
exact integers, and `IntegerPart` is `Listable`, threading over a vector of mixed
reals and rationals. Together with `FractionalPart` it splits any number as
`IntegerPart[x] + FractionalPart[x] == x`.
