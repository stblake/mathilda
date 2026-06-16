### Worked examples

```mathematica
In[1]:= IntegerLength[12345]
Out[1]= 5
```

```mathematica
In[1]:= IntegerLength[2^1000]
Out[1]= 302
```

```mathematica
In[1]:= IntegerLength[2^1000, 2]
Out[1]= 1001

In[2]:= IntegerLength[100!]
Out[2]= 158
```

### Notes

`IntegerLength[n]` returns the number of decimal digits of `n`, and
`IntegerLength[n, b]` the number of base-`b` digits. Working at arbitrary
precision, it answers questions that overflow fixed-width arithmetic: `2^1000`
has 302 decimal digits but exactly 1001 binary digits, and `100!` is a 158-digit
number. The sign of `n` is ignored and `IntegerLength[0]` is `0`.
