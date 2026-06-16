### Worked examples

```mathematica
In[1]:= Im[3 + 4 I]
Out[1]= 4

In[2]:= Im[7]
Out[2]= 0
```

`Im` resolves exact algebraic and transcendental values, including radicals of
negative numbers and branch-cut logarithms:

```mathematica
In[1]:= Im[Sqrt[-4]]
Out[1]= 2

In[2]:= Im[(1 + I)^10]
Out[2]= 32

In[3]:= Im[Log[-1]]
Out[3]= Pi
```

Being Listable, it threads over a list of numbers:

```mathematica
In[1]:= Im[{1 + 2 I, 3 - 4 I, 5}]
Out[1]= {2, -4, 0}
```

For values it cannot reduce in closed form, `Im` stays symbolic but still
yields to high-precision numerics — here the imaginary part of Γ(1 + i):

```mathematica
In[1]:= Im[Gamma[1 + I]]
Out[1]= Im[Gamma[1 + I]]

In[2]:= N[Im[Gamma[1 + I]], 40]
Out[2]= -0.15494982830181068512495513048388660519589
```

### Notes

`Im[z]` extracts the imaginary part of numeric `z`, giving 0 for real (or real-valued) arguments. It is Listable.
