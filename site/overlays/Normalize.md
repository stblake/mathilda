### Worked examples

```mathematica
In[1]:= Normalize[{3, 4}]
Out[1]= {3/5, 4/5}
```

Normalizing an integer vector whose norm is rational stays exact:

```mathematica
In[1]:= Normalize[{1, 2, 2}]
Out[1]= {1/3, 2/3, 2/3}
```

A symbolic vector normalizes to its general unit-vector formula:

```mathematica
In[1]:= Normalize[{a, b}]
Out[1]= {a/Sqrt[Abs[a]^2 + Abs[b]^2], b/Sqrt[Abs[a]^2 + Abs[b]^2]}
```

A complex scalar is divided by its modulus, giving a unit-modulus phase:

```mathematica
In[1]:= Normalize[3 + 4 I]
Out[1]= 3/5 + 4/5*I
```

`Normalize[expr, f]` normalizes with respect to any norm function — here `Total`
turns counts into a probability distribution:

```mathematica
In[1]:= Normalize[{1, 1}, Total]
Out[1]= {1/2, 1/2}
```

### Notes

`Normalize[v]` returns `v / Norm[v]`, the unit vector in the direction of `v`;
for a scalar (including a complex number) it returns `z / Abs[z]`. The two-argument
form `Normalize[expr, f]` divides by `f[expr]` instead, so any norm or aggregating
function may be supplied. Zero input is returned unchanged.
