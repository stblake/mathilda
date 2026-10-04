### Worked examples

```mathematica
In[1]:= InterpolatingPolynomial[{1, 4, 9}, x]  (* values at x = 1, 2, 3, in nested Newton form *)
```

```mathematica
In[1]:= Expand[InterpolatingPolynomial[{1, 8, 27, 64}, x]]  (* the cubes expand to x^3 *)
```

```mathematica
In[1]:= InterpolatingPolynomial[{{0, 0}, {1, 1}, {2, 4}}, x]  (* explicit {abscissa, value} pairs *)
```

```mathematica
In[1]:= InterpolatingPolynomial[{a, b, c}, x]  (* symbolic data gives a symbolic polynomial *)
```

```mathematica
In[1]:= InterpolatingPolynomial[{1, 4, 9, 16, 25}, x, Modulus -> 7]  (* the same fit over Z/7Z *)
```

### Notes

`InterpolatingPolynomial[{f1, f2, ...}, x]` gives the one polynomial reproducing
the values `fi` at `x = 1, 2, ...`, with `n` values yielding degree `n - 1`. The
default output is the nested **Newton–Horner** form — compact and numerically
stable — so `{1, 4, 9}` prints as `1 + (-1 + x) (1 + x)` rather than `x^2`; wrap in
`Expand` to see the monomial form.

Abscissae can be given explicitly as `{{xi, fi}, ...}`, and the values may be
exact, inexact, or symbolic — exact data give an exact polynomial. The
`{{{x1, y1, ...}, f1}, ...}` form with a variable list builds a multidimensional
interpolant of lowest total degree, and `Modulus -> n` performs the fit over the
integers mod `n`.
