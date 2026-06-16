### Worked examples

```mathematica
In[1]:= Re[3 + 4 I]
Out[1]= 3

In[2]:= Re[7]
Out[2]= 7
```

Being `Listable`, `Re` maps over a list of complex numbers, and it sees through
exact arithmetic — `(1 + I)^10 = 32 I` is purely imaginary, so its real part is
exactly `0`:

```mathematica
In[1]:= Re[{1 + I, 2 - 3 I, 5}]
Out[1]= {1, 2, 5}

In[2]:= Re[(1 + I)^10]
Out[2]= 0
```

It also rationalises quotients to extract an exact real part:

```mathematica
In[1]:= Re[1/(2 + 3 I)]
Out[1]= 2/13
```

### Notes

`Re[z]` extracts the real part of numeric `z`; a purely real argument is returned unchanged. It is Listable.
