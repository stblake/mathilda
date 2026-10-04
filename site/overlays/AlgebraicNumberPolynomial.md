### Worked examples

```mathematica
In[1]:= AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}], x]  (* c0 + c1 x from {1, 2} *)
```

```mathematica
In[1]:= AlgebraicNumberPolynomial[7, x]  (* a rational is the constant polynomial *)
```

### Notes

`AlgebraicNumberPolynomial[a, x]` gives the polynomial in `x` whose coefficients
are the stored coefficient list of the `AlgebraicNumber` object `a`: for
`a = AlgebraicNumber[theta, {c0, c1, …, cn}]` the result is
`c0 + c1 x + … + cn x^n`, and `a` is recovered by replacing `x` with `theta`. The
generator `theta` plays no part — this is a purely structural read of the
coefficient vector, with no field arithmetic.

An integer or rational `a` is the constant polynomial and is returned unchanged;
any other argument stays unevaluated. `Listable` and `Protected`.
