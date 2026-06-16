# Conjugate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Conjugate[z] gives the complex conjugate Re[z] - I Im[z] of numeric z; real and real-valued (Re/Im/Abs/Arg) arguments are returned unchanged.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_conjugate` folds the involution `Conjugate[Conjugate[z]] -> z` and treats the real-valued-by-construction heads `Re`, `Im`, `Abs`, `Arg` as fixed points. For a `Complex[re, im]` literal it returns `make_complex(re, -im)`; for real numerics (Integer/Real/Rational) and any expression that `is_numeric_real` (e.g. `Sqrt[2]`, `Pi`) it returns the argument unchanged. An expression that `complex_decompose` splits into concretely-numeric real/imag parts is conjugated as `re - im*I`. Symbolic inputs return `NULL` (a one-argument arity check emits `Conjugate::argx`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Conjugate[3 + 4 I]
Out[1]= 3 - 4*I

In[2]:= Conjugate[5]
Out[2]= 5
```

```mathematica
In[1]:= Conjugate[{1 + I, 2 - 3 I}]
Out[1]= {1 - I, 2 + 3*I}
```

```mathematica
In[1]:= Conjugate[(2 + I)/(1 - 3 I)]
Out[1]= -1/10 - 7/10*I
```

```mathematica
In[1]:= z Conjugate[z] /. z -> 3 + 4 I
Out[1]= 25
```

### Notes

`Conjugate[z]` returns `Re[z] - I Im[z]`; real arguments are returned unchanged. It is Listable.
