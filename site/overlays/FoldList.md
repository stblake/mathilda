### Worked examples

```mathematica
In[1]:= FoldList[Plus, 0, {1, 2, 3, 4}]
Out[1]= {0, 1, 3, 6, 10}
```

```mathematica
In[1]:= FoldList[Times, {1, 2, 3, 4, 5}]
Out[1]= {1, 2, 6, 24, 120}
```

```mathematica
In[1]:= FoldList[Max, {3, 1, 4, 1, 5, 9, 2, 6}]
Out[1]= {3, 3, 4, 4, 5, 9, 9, 9}
```

```mathematica
In[1]:= FoldList[(#1 + #2)/2 &, 0, {1, 1, 1, 1}]
Out[1]= {0, 1/2, 3/4, 7/8, 15/16}
```

### Notes

`FoldList[f, x, list]` returns every intermediate accumulator, giving a
length-`n+1` list. `FoldList[Plus, 0, l]` produces cumulative sums and
`FoldList[Times, l]` the running partial factorials. `FoldList[Max, l]`
yields the running maximum (a prefix-scan idiom), and the exact-arithmetic
midpoint iteration `(#1 + #2)/2 &` converges on the dyadic rationals
`1 - 2^-k`, kept exact as `Rational`s rather than floats.
