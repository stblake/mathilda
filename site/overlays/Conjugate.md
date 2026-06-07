### Worked examples

```mathematica
In[1]:= Conjugate[3 + 4 I]
Out[1]= 3 - 4*I

In[2]:= Conjugate[5]
Out[2]= 5

In[3]:= Conjugate[{1 + I, 2 - 3 I}]
Out[3]= {1 - I, 2 + 3*I}
```

### Notes

`Conjugate[z]` returns `Re[z] - I Im[z]`; real arguments are returned unchanged. It is Listable.
