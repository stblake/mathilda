### Worked examples

```mathematica
In[1]:= ProductLog[0]
```

```mathematica
In[1]:= ProductLog[E]  (* since e = 1 * e^1, the principal value is 1 *)
```

```mathematica
In[1]:= ProductLog[-1/E]  (* the branch point where the k = 0 and k = -1 branches meet *)
```

```mathematica
In[1]:= ProductLog[-Pi/2]
```

```mathematica
In[1]:= N[ProductLog[1], 30]  (* the omega constant *)
```

```mathematica
In[1]:= ProductLog[-1, -0.2]  (* the k = -1 branch *)
```

```mathematica
In[1]:= D[ProductLog[x], x]
```

### Notes

`ProductLog[z]` is the Lambert W function: the principal branch `W_0(z)`, the
solution `w` of `z = w e^w`. `ProductLog[k, z]` selects the `k`-th branch (`k`
an explicit integer); `k = 0` is the principal branch and `k = -1` is the other
real branch.

A handful of exact values are returned directly (`0`, `E`, `-1/E`, `-Pi/2`);
everything else is evaluated numerically by a Halley iteration in
arbitrary-precision complex arithmetic, returning a real result on the
real-valued domains (`W_0` for `x >= -1/e`, `W_{-1}` for `-1/e <= x < 0`) and a
`Complex[..]` otherwise. `N[ProductLog[1], 30]` is the omega constant `Omega`,
the solution of `Omega e^Omega = 1`. The derivative is
`ProductLog[z]/(z (1 + ProductLog[z]))`.
