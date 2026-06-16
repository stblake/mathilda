# Im

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Im[z] gives the imaginary part of numeric z, and 0 for real or real-valued (Re/Im/Abs/Arg) arguments.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_im` returns the imaginary part. It returns `0` for the real-valued-by-construction heads (`Re`/`Im`/`Abs`/`Arg`) and for any real numeric kind (Integer/Real/Rational/MPFR), copies the second component of a `Complex[re, im]` literal, and for a general expression runs `complex_decompose` (a recursive Plus/Times walk that propagates `Complex` literals through complex multiplication) — returning the imaginary part only when both decomposed parts are concretely numeric (`is_numeric_real`). Otherwise `NULL`, leaving the symbolic head in place. `Re`/`ReIm` in the same file share this machinery.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Im[3 + 4 I]
Out[1]= 4

In[2]:= Im[7]
Out[2]= 0
```

`Im` resolves exact algebraic and transcendental values, including radicals of
negative numbers and branch-cut logarithms:

```mathematica
In[1]:= Im[Sqrt[-4]]
Out[1]= 2

In[2]:= Im[(1 + I)^10]
Out[2]= 32

In[3]:= Im[Log[-1]]
Out[3]= Pi
```

Being Listable, it threads over a list of numbers:

```mathematica
In[1]:= Im[{1 + 2 I, 3 - 4 I, 5}]
Out[1]= {2, -4, 0}
```

For values it cannot reduce in closed form, `Im` stays symbolic but still
yields to high-precision numerics — here the imaginary part of Γ(1 + i):

```mathematica
In[1]:= Im[Gamma[1 + I]]
Out[1]= Im[Gamma[1 + I]]

In[2]:= N[Im[Gamma[1 + I]], 40]
Out[2]= -0.15494982830181068512495513048388660519589
```

### Notes

`Im[z]` extracts the imaginary part of numeric `z`, giving 0 for real (or real-valued) arguments. It is Listable.
