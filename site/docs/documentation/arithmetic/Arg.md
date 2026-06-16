# Arg

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Arg[z] gives the argument (phase angle in (-Pi, Pi]) of numeric z; 0 for nonnegative reals, Pi for negative reals.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_arg` returns the phase angle in `(-Pi, Pi]`. A pure MPFR real folds to exact `0` or `Pi` by sign. For a `Complex[re, im]` whose parts are exact (Integer/Rational), it recognises the special directions and returns exact multiples of `Pi`: `0` for positive reals, `Pi` for negatives, `±Pi/2` on the imaginary axis, and `±Pi/4`, `±3Pi/4` on the diagonals; otherwise it returns the symbolic `ArcTan[re, im]`. When either component carries MPFR it evaluates `mpfr_atan2` at the combined precision; an inexact machine `Real` falls through to the libm `atan2(im, re)`. Symbolic inputs return `NULL`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/complex.c`](https://github.com/stblake/mathilda/blob/main/src/complex.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Arg[1]
Out[1]= 0

In[2]:= Arg[-1]
Out[2]= Pi

In[3]:= Arg[1 + I]
Out[3]= 1/4 Pi
```

Powers wrap the phase into the principal branch:

```mathematica
In[1]:= Arg[(1 + I)^10]
Out[1]= 1/2 Pi

In[2]:= Arg[-2 + 2 I]
Out[2]= 3/4 Pi
```

Stepping the argument of `(1 + I)^k` traces the unit circle and folds back at `Pi`:

```mathematica
In[1]:= Table[Arg[(1 + I)^k], {k, 0, 8}]
Out[1]= {0, 1/4 Pi, 1/2 Pi, 3/4 Pi, Pi, -3/4 Pi, -1/2 Pi, -1/4 Pi, 0}
```

Generic complex numbers evaluate to high precision:

```mathematica
In[1]:= N[Arg[2 + 3 I], 40]
Out[1]= 0.98279372324732906798571061101466601449686
```

### Notes

`Arg[z]` gives the phase angle in the range `(-Pi, Pi]`: 0 for positive reals, `Pi` for negative reals.
