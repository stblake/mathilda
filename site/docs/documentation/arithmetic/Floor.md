# Floor

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Floor[x]
    gives the greatest integer less than or equal to x.
Floor[x, a]
    gives the greatest multiple of a less than or equal to x.
Floor is Listable. Exact (Integer / BigInt / Rational) inputs return
exact integers; Real / MPFR inputs are rounded toward -Infinity at
the input precision; symbolic inputs stay unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_floor` calls the shared `do_piecewise(res, OP_FLOOR, "Floor", true)` dispatcher (the `true` enables the two-argument `Floor[x, a]` → `a Floor[x/a]` unit-quantization form). The per-element kernel `do_piecewise_1` handles each numeric type exactly: `EXPR_INTEGER`/`EXPR_BIGINT` are already integers and pass through; rationals are floored exactly via GMP; `EXPR_REAL` uses C `floor()` cast to `int64_t`; `EXPR_MPFR` uses `mpfr_floor` then `mpfr_get_z` into an `mpz_t` (normalized to int/bigint) so arbitrarily large values never silently truncate; `±Infinity` pass through. `Floor` is registered `LISTABLE | NUMERICFUNCTION | PROTECTED`, so threading over lists is handled by the generic evaluator before the builtin runs. Non-numeric arguments return NULL (left symbolic).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/piecewise.c`](https://github.com/stblake/mathilda/blob/main/src/piecewise.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Floor[7/2]
Out[1]= 3
```

```mathematica
In[1]:= Floor[-2.3]
Out[1]= -3
```

```mathematica
In[1]:= Floor[17, 5]
Out[1]= 15
```

```mathematica
In[1]:= Floor[{2.7, -2.7, 5/2, 11/3}]
Out[1]= {2, -3, 2, 3}
```

```mathematica
In[1]:= Floor[N[Pi, 40] 10^20]
Out[1]= 314159265358979323846
```

### Notes

`Floor[x]` rounds toward `-Infinity`; the two-argument `Floor[x, a]` gives the
greatest multiple of `a` not exceeding `x`. Exact inputs return exact
integers, and `Floor` is `Listable` so it threads over a vector of mixed
reals and rationals. The last example extracts the first 21 significant
digits of `π` exactly by flooring a 40-digit MPFR value scaled by `10^20`.
