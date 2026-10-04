### Worked examples

```mathematica
In[1]:= HammingDistance[{1, 2, 3, 4}, {1, 0, 3, 0}]  (* two positions differ *)
```

```mathematica
In[1]:= HammingDistance[{1, 0, 1, 1}, {1, 1, 0, 1}]  (* bit-vector distance *)
```

```mathematica
In[1]:= HammingDistance["2718281828", "3141592653"]  (* compared character by character *)
```

### Notes

`HammingDistance[a, b]` counts the positions at which two equal-length sequences
differ. Both arguments must be strings, or both lists; a length mismatch leaves
the call unevaluated (unlike `EditDistance`, which handles differing lengths).
Comparison is by equality, so it works on bit vectors, digit strings, or lists of
arbitrary expressions. Strings are compared byte by byte, so a multi-byte UTF-8
character spans several positions.
