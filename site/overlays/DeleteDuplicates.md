### Worked examples

```mathematica
In[1]:= DeleteDuplicates[{1,2,1,3,2,4}]
Out[1]= {1, 2, 3, 4}
```

```mathematica
In[1]:= DeleteDuplicates[{a,b,a,c,b,d,a}]
Out[1]= {a, b, c, d}
```

```mathematica
In[1]:= DeleteDuplicates[{1,2,3,4,5,6}, Mod[#1,3]==Mod[#2,3]&]
Out[1]= {1, 2, 3}
```

```mathematica
In[1]:= DeleteDuplicates[{1,-1,2,-2,3,-3}, Abs[#1]==Abs[#2]&]
Out[1]= {1, 2, 3}
```

### Notes

`DeleteDuplicates[list]` keeps the first occurrence of each distinct element and preserves the original order. With a two-argument equivalence test `DeleteDuplicates[list, test]` two elements are treated as duplicates when `test[a, b]` is `True`, so you can deduplicate by an arbitrary relation rather than literal equality — picking one representative per residue class modulo 3, or one per absolute value, in the examples above.
