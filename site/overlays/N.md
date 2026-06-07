---
status: Stable
references:
---
### Worked examples

```mathematica
In[1]:= N[Pi, 20]
Out[1]= 3.14159265358979323846
```

```mathematica
In[1]:= N[Sqrt[2]]
Out[1]= 1.41421
```

```mathematica
In[1]:= N[2/7, 15]
Out[1]= 0.2857142857142856
```

```mathematica
In[1]:= N[E]
Out[1]= 2.71828
```

### Notes

`N[expr]` gives a machine-precision floating-point value, displayed to about six
significant digits. `N[expr, d]` requests approximately `d` digits of precision,
computed via arbitrary-precision arithmetic (so `N[Pi, 20]` returns the constant
to 20 digits). Exact inputs such as `Sqrt[2]`, `Pi`, `E`, and rationals are
converted to their numeric approximations. Note that machine-precision results
print at the default short width even when more digits are internally available.
