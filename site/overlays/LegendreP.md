### Worked examples

```mathematica
In[1]:= LegendreP[2, x]
```

```mathematica
In[1]:= Table[LegendreP[n, x], {n, 0, 4}]
```

```mathematica
In[1]:= LegendreP[4, 0]  (* an even-degree value at the origin *)
```

```mathematica
In[1]:= LegendreP[n, 1]  (* every P_n equals 1 at the endpoint *)
```

```mathematica
In[1]:= LegendreP[2, 1, x]  (* the associated function P_2^1 *)
```

```mathematica
In[1]:= Integrate[LegendreP[2, x] LegendreP[3, x], {x, -1, 1}]  (* orthogonality on the interval *)
```

```mathematica
In[1]:= N[LegendreP[1.5, 0.3]]  (* non-integer order via the Gauss 2F1 series *)
```

```mathematica
In[1]:= D[LegendreP[3, x], x]
```

### Notes

For an exact integer order `LegendreP[n, x]` is the explicit Legendre polynomial
`P_n(x)`, built from the three-term recurrence with exact rational coefficients
(orders up to `n = 2000`); `LegendreP[n, 1]` is `1` for every order. The
Legendre polynomials are orthogonal on `[-1, 1]`, so the integral of a product
of two of different degree vanishes.

A non-integer order is evaluated numerically (when an argument is inexact)
through the Gauss hypergeometric series `P_n(x) = 2F1(-n, n+1; 1; (1-x)/2)`,
valid for `|(1-x)/2| < 1`.

`LegendreP[n, m, x]` is the associated Legendre function, and
`LegendreP[n, m, a, x]` selects the Legendre function of type `a` (one of 1, 2,
3). The derivatives of the integer-order polynomials follow from the explicit
form.
