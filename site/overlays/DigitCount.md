### Worked examples

```mathematica
In[1]:= DigitCount[1234567890]
Out[1]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
```

```mathematica
In[1]:= DigitCount[100!]
Out[1]= {15, 19, 10, 10, 14, 19, 7, 14, 20, 30}
```

```mathematica
In[1]:= DigitCount[255, 2]
Out[1]= {8, 0}
```

```mathematica
In[1]:= DigitCount[2^1000, 10, 0]
Out[1]= 28
```

### Notes

`DigitCount[n]` tallies how often each digit appears in the base-10 expansion of
`n`, in the order `1, 2, ..., 9, 0` (note that `0` is reported last). The pandigital
number `1234567890` therefore gives ten ones. The distribution of the 158 digits of
`100!` shows the histogram for a large factorial. With a base argument,
`DigitCount[n, b]` works in any radix — `DigitCount[255, 2] = {8, 0}` recovers the
binary population count (255 is `11111111`). The three-argument form
`DigitCount[n, b, d]` returns just the count of digit `d`; here `2^1000` (a 302-digit
number) contains 28 zeros. The sign of `n` is ignored and `DigitCount[0]` is all
zeros.
