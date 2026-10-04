### Worked examples

```mathematica
In[1]:= QPochhammer[a, q, 3]  (* the finite q-shifted factorial *)
```

```mathematica
In[1]:= QPochhammer[a, q, 0]
```

```mathematica
In[1]:= QPochhammer[2, 1/3, 4]  (* exact rational arguments collapse to a rational *)
```

```mathematica
In[1]:= QPochhammer[0.5, 0.5]  (* the infinite form at machine precision *)
```

### Notes

`QPochhammer[a, q, n]` is the finite q-Pochhammer symbol (q-shifted factorial)
`prod_{k=0}^{n-1} (1 - a q^k)`, with `QPochhammer[a, q, 0] = 1`. For a
non-negative integer `n` the product is formed and evaluated — exactly for exact
`a, q`, numerically when they are inexact. A symbolic or non-integer `n` is left
unevaluated, which is what `Product` relies on to return the symbol as a closed
form.

The two-argument `QPochhammer[a, q]` is the infinite product `(a; q)_inf`,
evaluated for machine-real `a, q` with `|q| < 1`. That infinite form is a
machine-precision kernel only, so a request for extra digits via `N[..., p]`
still returns a machine-precision value.
