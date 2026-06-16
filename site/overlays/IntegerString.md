### Worked examples

```mathematica
In[1]:= IntegerString[255, 16]
Out[1]= "ff"
```

```mathematica
In[1]:= IntegerString[255, 2]
Out[1]= "11111111"

In[2]:= IntegerString[42, 2, 16]
Out[2]= "0000000000101010"
```

```mathematica
In[1]:= IntegerString[123456789, 36]
Out[1]= "21i3v9"
```

```mathematica
In[1]:= IntegerString[3^50, 16]
Out[1]= "980553f0db2fd09de3c9"
```

### Notes

`IntegerString[n, b]` renders `n` as a base-`b` string, using the digits `0-9`
and then the letters `a-z` for values 10 through 35 (so the maximum base is 36).
A third argument left-pads with zeros to a fixed width. Because the conversion
runs at arbitrary precision, large numbers like `3^50` are formatted exactly in
hexadecimal. The sign of `n` is discarded.
