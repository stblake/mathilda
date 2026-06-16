### Worked examples

```mathematica
In[1]:= Norm[{3, 4}]
Out[1]= 5
```

A symbolic 2-norm is kept exact, written with `Abs` so it is valid for complex
components too:

```mathematica
In[1]:= Norm[{a, b, c}]
Out[1]= Sqrt[Abs[a]^2 + Abs[b]^2 + Abs[c]^2]
```

The same vector under different p-norms — the 1-norm, the max (Infinity) norm,
and the exact 3-norm:

```mathematica
In[1]:= Norm[{1, 2, 3, 4}, 1]
Out[1]= 10

In[2]:= Norm[{1, 2, 3, 4}, Infinity]
Out[2]= 4

In[3]:= Norm[{1, 2, 3, 4}, 3]
Out[3]= 10^(2/3)
```

The norm of a complex vector, taken to 40 significant digits:

```mathematica
In[1]:= N[Norm[{1, 2, 3}], 40]
Out[1]= 3.7416573867739413855837487323165493017559
```

### Notes

`Norm[expr]` gives the 2-norm of a scalar or vector; `Norm[expr, p]` gives the
p-norm (`Abs` for scalars; `(Sum |xi|^p)^(1/p)` for vectors with `1 <= p <
Infinity`; `Max[Abs[expr]]` for `p == Infinity`). Results stay exact for exact
input, and symbolic components are written with `Abs` so the formula is valid for
complex entries.
