# ReIm

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ReIm[z] gives {Re[z], Im[z]}, the real and imaginary parts of numeric z as a list; real-valued arguments give {z, 0}.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_reim` returns `{Re[z], Im[z]}` as a two-element `List`. It mirrors `builtin_re`/`builtin_im`: `{f[z], 0}` for real-valued head calls (`Re`/`Im`/`Abs`/`Arg`), `{re, im}` for a `Complex[re, im]` literal, `{x, 0}` for real numeric kinds (Integer/Real/Rational), and `{re, im}` when `complex_decompose` yields numeric parts. Returns `NULL` (unevaluated) for genuinely symbolic input.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ReIm[3 + 4 I]
Out[1]= {3, 4}
```

It splits exact powers and quotients into their real/imaginary components —
`(2 + I)^3 = 2 + 11 I`, and a complex division reduces to integers:

```mathematica
In[1]:= ReIm[(2 + I)^3]
Out[1]= {2, 11}

In[2]:= ReIm[(3 + 4 I)/(1 - 2 I)]
Out[2]= {-1, 2}
```

On a transcendental argument it returns the numeric pair — here Euler's formula
`e^(iπ/4)` to 20 digits, the real and imaginary parts each `1/√2`:

```mathematica
In[1]:= ReIm[N[E^(I Pi/4), 20]]
Out[1]= {0.707106781186547524409, 0.707106781186547524395}
```

### Notes

`ReIm[z]` is shorthand for `{Re[z], Im[z]}`; a real-valued argument gives `{z, 0}`.
