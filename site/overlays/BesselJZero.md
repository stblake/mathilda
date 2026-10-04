### Worked examples

```mathematica
In[1]:= BesselJZero[0, 1]  (* no numeric path yet -- stays symbolic *)
```

```mathematica
In[1]:= Head[BesselJZero[n, k]]  (* a bare, unevaluated head *)
```

```mathematica
In[1]:= Product[1 - x^2/BesselJZero[n, k]^2, {k, 1, Infinity}]  (* the Hadamard product it exists to name *)
```

### Notes

`BesselJZero[n, k]` names the `k`-th positive zero of `BesselJ[n, x]`. It is currently
a **symbolic placeholder**: it carries no numeric evaluator, so every call — even
`BesselJZero[0, 1]` — stays unevaluated.

Its reason for existing today is the infinite-product recogniser in the `Product`
subsystem, which matches the canonical Hadamard product

```
Product[1 - x^2/BesselJZero[n, k]^2, {k, 1, Infinity}] = Gamma[1 + n] (2/x)^n BesselJ[n, x]
```

and so reproduces the closed Bessel form directly. `BesselJZero` is `Listable`,
`NumericFunction` and `Protected`.
