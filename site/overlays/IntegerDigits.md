### Worked examples

```mathematica
In[1]:= IntegerDigits[12345]
Out[1]= {1, 2, 3, 4, 5}
```

```mathematica
In[1]:= IntegerDigits[255, 16]
Out[1]= {15, 15}

In[2]:= IntegerDigits[255, 2]
Out[2]= {1, 1, 1, 1, 1, 1, 1, 1}

In[3]:= IntegerDigits[5, 2, 8]
Out[3]= {0, 0, 0, 0, 0, 1, 0, 1}
```

```mathematica
In[1]:= Total[IntegerDigits[2^100]]
Out[1]= 115
```

```mathematica
In[1]:= IntegerDigits[100!, 10][[1 ;; 5]]
Out[1]= {9, 3, 3, 2, 6}
```

### Notes

`IntegerDigits[n]` returns the decimal digits of `n` most-significant first;
`IntegerDigits[n, b]` works in base `b`, and a third argument left-pads (or
truncates to the least-significant) to a fixed length. Because Mathilda carries
arbitrary-precision integers, the digit list of a giant number such as `2^100`
or `100!` is exact — handy for digit-sum problems and divisibility tricks. The
sign of `n` is ignored.
