# ArcTan

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ArcTan[z]
    gives the principal inverse tangent of z, in (-Pi/2, Pi/2).
ArcTan[y, x]
    gives the argument of the complex number x + I y, in (-Pi, Pi]
    (two-argument atan2 form).
ArcTan is Listable.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_arctan` (`src/trig.c`) handles one- and two-argument forms. For `ArcTan[z]`: (1) `odd_fold` for oddness; (2) `trig_i_fold` for `ArcTan[I y] -> I ArcTanh[y]`; (3) exact inversion via `exact_arctan`, scanning n in `[-d/2, d/2]` for d in {1,2,3,4,5,6,10,12} against `exact_tan(n,d)` and returning `n/d * Pi`; (4) numeric fallback MPFR `mpfr_atan`/`mpfr_complex_atan`, else `get_approx` + C99 `catan`. For the two-argument `ArcTan[x, y]` (argument of x + I y): integer inputs are resolved exactly by quadrant (the axes and the four diagonals `±1` map to `0, ±Pi/2, Pi, ±Pi/4, ±3Pi/4`); real inputs use MPFR `mpfr_atan2` or C99 `atan2`. Indeterminate `ArcTan[0,0]` and unhandled cases return `NULL`.

**Data structures.** `Expr*` trees; the single-arg exact path reuses the forward `exact_tan` table.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/trig.c`](https://github.com/stblake/mathilda/blob/main/src/trig.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= ArcTan[1]
Out[1]= 1/4 Pi

In[2]:= ArcTan[Sqrt[3]]
Out[2]= 1/3 Pi
```

```mathematica
In[1]:= ArcTan[-1, 1]
Out[1]= 3/4 Pi

In[2]:= ArcTan[{0, 1}]
Out[2]= {0, 1/4 Pi}
```

A pure-imaginary argument crosses over to the inverse hyperbolic tangent:

```mathematica
In[1]:= ArcTan[I/2]
Out[1]= I ArcTanh[1/2]
```

The Taylor series recovers the Gregory series, and the antiderivative is closed-form:

```mathematica
In[1]:= Series[ArcTan[x], {x, 0, 9}]
Out[1]= x - 1/3 x^3 + 1/5 x^5 - 1/7 x^7 + 1/9 x^9 + O[x]^10

In[2]:= Integrate[ArcTan[x], x]
Out[2]= 1/2 (2 x ArcTan[x] - Log[1 + x^2])
```

Machin's 1706 formula then evaluates Pi to 40 digits:

```mathematica
In[1]:= N[16 ArcTan[1/5] - 4 ArcTan[1/239], 40]
Out[1]= 3.1415926535897932384626433832795028841975
```

### Notes

`ArcTan[z]` gives the principal inverse tangent, in `(-Pi/2, Pi/2)`. The two-argument form `ArcTan[x, y]` is the quadrant-aware `atan2`, giving the argument of `x + I y` in `(-Pi, Pi]`. `ArcTan` is Listable.
