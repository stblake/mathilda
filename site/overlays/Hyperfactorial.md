### Worked examples

```mathematica
In[1]:= Hyperfactorial[4]  (* the product 1^1 2^2 3^3 4^4 *)
```

```mathematica
In[1]:= Hyperfactorial[0]
```

```mathematica
In[1]:= Table[Hyperfactorial[n], {n, 0, 5}]
```

```mathematica
In[1]:= Hyperfactorial[1/2]  (* exact non-integer orders stay symbolic *)
```

```mathematica
In[1]:= N[Hyperfactorial[5/2], 20]  (* the Barnes-G continuation, to 20 digits *)
```

### Notes

`Hyperfactorial[n]` is the product `prod_{k=1}^{n} k^k`, with
`Hyperfactorial[0] = Hyperfactorial[1] = 1`. A non-negative integer order is
computed exactly in arbitrary-precision integer arithmetic (capped at
`n = 20000`), so the values grow very fast: `1, 1, 4, 108, 27648, 86400000, ...`.

A non-integer or complex order is evaluated only numerically (under `N`),
through the Barnes-G analytic continuation `Hyperfactorial[z] = Gamma[z+1]^z /
BarnesG[z+1]`, which reproduces the integer values and extends the function to
the whole plane. Exact non-integer and negative-integer orders are left
symbolic.

The head exists partly so that `Product` can recognise `prod k^k` and return it
in closed form.
