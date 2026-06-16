### Worked examples

```mathematica
In[1]:= x = 10
Out[1]= 10

In[2]:= Block[{x = 3}, x^2]
Out[2]= 9

In[3]:= x
Out[3]= 10
```

```mathematica
In[1]:= Block[{n = 5}, Sum[k, {k, 1, n}]]
Out[1]= 15
```

### Notes

`Block` localizes the *values* of the listed symbols: it temporarily resets them for the duration of the body and restores the outer values on exit. Unlike `Module`, it does not rename the symbols, so the same global `x` is reused with a saved-and-restored value.
