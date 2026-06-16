### Worked examples

```mathematica
In[1]:= FromDigits[{1, 2, 3, 4}]
Out[1]= 1234
```

```mathematica
In[1]:= FromDigits[{1, 0, 1, 1}, 2]
Out[1]= 11
```

```mathematica
In[1]:= FromDigits["deadbeef", 16]
Out[1]= 3735928559
```

```mathematica
In[1]:= FromDigits[IntegerDigits[2^100], 10]
Out[1]= 1267650600228229401496703205376
```

```mathematica
In[1]:= FromDigits[{d2, d1, d0}, b]
Out[1]= d0 + b d1 + b^2 d2
```

### Notes

`FromDigits` evaluates a digit list most-significant-first via Horner's method.
A second argument fixes the base: `{1, 0, 1, 1}` in base 2 is `11`, and the hex
string `"deadbeef"` (letters `a`-`z` denoting 10-35) is `3735928559`. Because the
arithmetic is exact arbitrary-precision, round-tripping a 31-digit number through
`IntegerDigits`/`FromDigits` reproduces it exactly. With symbolic digits or base,
the Horner recurrence is returned as a polynomial, `d0 + b d1 + b^2 d2` — the
inverse construction to `IntegerDigits`.
