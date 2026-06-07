# Sqrt

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Sqrt[z]
    represents the principal square root of z.
Sqrt is Listable. Sqrt[z] is canonicalised to Power[z, 1/2]; perfect
integer / rational squares reduce to exact form, negative real inputs
yield I * Sqrt[-x], and numeric inputs (Real / MPFR / Complex) are
evaluated directly. Branch cut along the negative real axis.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_sqrt` is a thin wrapper: it rewrites `Sqrt[x]` to `Power[x, Rational[1, 2]]` (via `make_rational(1, 2)`) and returns that, letting the full `Power` machinery handle all simplification (exact perfect squares, `Sqrt[8] -> 2 Sqrt[2]` radical extraction, numeric/MPFR evaluation, infinity algebra). `Sqrt` carries `LISTABLE | NUMERICFUNCTION | PROTECTED`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on square-free factorization and radical simplification.
- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on integer square roots.
- Source: [`src/power.c`](https://github.com/stblake/mathilda/blob/main/src/power.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Sqrt[50]
Out[1]= 5 Sqrt[2]
```

```mathematica
In[1]:= Sqrt[18]
Out[1]= 3 Sqrt[2]
```

```mathematica
In[1]:= Sqrt[-9]
Out[1]= 3*I
```

```mathematica
In[1]:= Sqrt[1/4]
Out[1]= 1/2
```

### Notes

`Sqrt[n]` is `Power[n, 1/2]`, so it inherits the perfect-square extraction logic:
the largest square factor is pulled out of the radical, reducing `Sqrt[50]` to
`5 Sqrt[2]` and `Sqrt[18]` to `3 Sqrt[2]`. Perfect squares collapse fully, and a
perfect-square rational like `1/4` yields the exact rational `1/2`. Negative
arguments produce a pure imaginary result, with `Sqrt[-9]` giving `3*I`. The
surd is kept in symbolic form when no square factor remains, e.g. `Sqrt[2]`.
