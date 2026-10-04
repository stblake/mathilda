### Worked examples

```mathematica
In[1]:= ListGradient[{1, 2, 4, 7, 11}]  (* central differences inside, one-sided at the ends *)
```

```mathematica
In[1]:= ListGradient[{a, b, c, d}]  (* symbolic input gives a symbolic gradient *)
```

```mathematica
In[1]:= ListGradient[{1., 2., 4., 7., 11.}]  (* machine input rides the float buffer *)
```

```mathematica
In[1]:= ListGradient[{1, 2, 4, 7, 11}, 2]  (* uniform spacing h = 2 *)
```

```mathematica
In[1]:= ListGradient[{0, 1, 4, 9, 16}, {0, 1, 2, 4, 8}]  (* a non-uniform coordinate grid *)
```

```mathematica
In[1]:= ListGradient[{{1, 2, 6}, {3, 4, 5}}]  (* rank 2: one array per axis *)
```

```mathematica
In[1]:= ListGradient[{1, 2, 4, 7, 11}, Method -> "Forward"]  (* forward differences everywhere *)
```

### Notes

`ListGradient[f]` is Mathilda's port of `numpy.gradient`: it estimates the
derivative of a sampled array by finite differences. By default it uses a
second-order central difference in the interior and a first-order one-sided
difference at the two extreme endpoints, so the result has the same shape as `f`.
A vector returns a vector; a rank-*k* array returns `{g_1, ..., g_k}`, one array
per axis.

Spacing is given by the optional second argument: a scalar `h` is uniform on
every axis, `{s_1, ..., s_k}` gives one spec per axis, and a coordinate list the
same length as the axis samples a **non-uniform** grid — handled with no special
case, because the underlying Fornberg weights are exact on arbitrarily spaced
nodes. Options `Method`, `DifferenceOrder`, `WindowLength` and `Axis` select the
one-sided variants, the accuracy order, the stencil size, and which axes to
differentiate.

Exactness follows the input: an integer or rational grid yields exact Rationals
(`(f[i+1]-f[i-1])/2` and friends), a symbolic grid yields symbolic weights, and a
machine-real (or complex) array runs straight on the packed float buffer and
stays packed. Because an integer gradient is Rational and no int64 buffer holds
one, an exact integer input is routed to the List path rather than the buffer.
The rank-1 real-array form also lowers inside `Compile[]`. Each differentiated
axis must have length at least 2.
