### Worked examples

```mathematica
In[1]:= Log[E]
Out[1]= 1

In[2]:= Log[E^2]
Out[2]= 2

In[3]:= Log[2, 8]
Out[3]= 3

In[4]:= Log[{1, E, E^2}]
Out[4]= {0, 1, 2}
```

```mathematica
In[1]:= Log[-1]
Out[1]= I Pi
```

```mathematica
In[1]:= D[Log[Sin[x]], x]
Out[1]= Cot[x]
```

```mathematica
In[1]:= Series[Log[1 + x], {x, 0, 6}]
Out[1]= x - 1/2 x^2 + 1/3 x^3 - 1/4 x^4 + 1/5 x^5 - 1/6 x^6 + O[x]^7
```

```mathematica
In[1]:= N[Log[2], 40]
Out[1]= 0.69314718055994530941723212145817656807549
```

```mathematica
In[1]:= N[Log[1 + I], 40]
Out[1]= 0.34657359027997265470861606072908828403779 + 0.78539816339744830961566084581987572104928*I
```

### Notes

`Log[z]` is the principal natural logarithm; `Log[b, z]` gives the base-`b` logarithm `Log[z]/Log[b]`. Log is Listable.
