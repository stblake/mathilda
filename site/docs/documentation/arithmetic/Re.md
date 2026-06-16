# Re

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Re[z] gives the real part of numeric z; Re[Re[z]], Re[Im[z]], Re[Abs[z]], Re[Arg[z]] fold since those heads are real-valued.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_re` returns the real part. It returns the argument itself for real numeric kinds (`EXPR_INTEGER`/`EXPR_REAL`/`EXPR_MPFR`/Rational) and for real-valued head calls (`Re`/`Im`/`Abs`/`Arg`, via `is_real_valued_head_call`); for a `Complex[re, im]` literal it returns `re`; and for an expression `complex_decompose` splits into numeric real/imaginary parts, it returns the real part. Otherwise (genuinely symbolic) it returns `NULL` and the call stays unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Re[3 + 4 I]
Out[1]= 3

In[2]:= Re[7]
Out[2]= 7
```

Being `Listable`, `Re` maps over a list of complex numbers, and it sees through
exact arithmetic — `(1 + I)^10 = 32 I` is purely imaginary, so its real part is
exactly `0`:

```mathematica
In[1]:= Re[{1 + I, 2 - 3 I, 5}]
Out[1]= {1, 2, 5}

In[2]:= Re[(1 + I)^10]
Out[2]= 0
```

It also rationalises quotients to extract an exact real part:

```mathematica
In[1]:= Re[1/(2 + 3 I)]
Out[1]= 2/13
```

### Notes

`Re[z]` extracts the real part of numeric `z`; a purely real argument is returned unchanged. It is Listable.
