### Worked examples

```mathematica
In[1]:= Solve[x^2 + 1 == 0, x, Modulus -> 5]  (* roots in Z/5Z *)
```

```mathematica
In[1]:= Factor[x^2 + 1, Modulus -> 5]  (* reducible over GF(5): (x+2)(x+3) *)
```

```mathematica
In[1]:= Solve[x^2 == 2, x, Modulus -> 7]  (* 2 is a quadratic residue mod 7 *)
```

### Notes

`Modulus -> p` is an **option** that moves an operation into the finite ring `Z/pZ`
(or the field `GF(p)` for prime `p`). It is read by `Solve`, `Factor`, `PolynomialGCD`,
`PolynomialReduce`, `GroebnerBasis`, `Reduce` and related polynomial heads; `Modulus`
itself is an inert, `Protected` symbol with no value of its own.

For `Solve`, a single-variable polynomial equation is solved by residue enumeration,
returning `{{x -> r}, ...}` with `r` ascending in `[0, p)`, for `2 <= p <= 100000`;
systems, multivariable specs, non-polynomial equations, or an out-of-range modulus
leave `Solve` unevaluated. The polynomial heads generally require `p` prime (a composite
is reported unsupported rather than silently mis-factored); the default `Modulus -> 0`
is the ordinary integer ring.
