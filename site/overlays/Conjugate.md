### Worked examples

```mathematica
In[1]:= Conjugate[3 + 4 I]
Out[1]= 3 - 4*I

In[2]:= Conjugate[5]
Out[2]= 5
```

```mathematica
In[1]:= Conjugate[{1 + I, 2 - 3 I}]
Out[1]= {1 - I, 2 + 3*I}
```

```mathematica
In[1]:= Conjugate[(2 + I)/(1 - 3 I)]
Out[1]= -1/10 - 7/10*I
```

```mathematica
In[1]:= z Conjugate[z] /. z -> 3 + 4 I
Out[1]= 25
```

### Notes

`Conjugate[z]` returns `Re[z] - I Im[z]`; real arguments are returned unchanged. It is Listable.
