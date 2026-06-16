### Worked examples

```mathematica
In[1]:= Count[{1,2,1,3,1},1]
Out[1]= 3
```

```mathematica
In[1]:= Count[{1,2,3,4,5,6},_?EvenQ]
Out[1]= 3
```

```mathematica
In[1]:= Count[Range[100], _?PrimeQ]
Out[1]= 25
```

```mathematica
In[1]:= Count[{a, b, {c, a}, a, {a, {a}}}, a, Infinity]
Out[1]= 5
```

```mathematica
In[1]:= Count[IntegerDigits[2^100], _?(# > 5 &)]
Out[1]= 11
```

### Notes

The second argument is a pattern, so `Count[list, _?EvenQ]` counts elements satisfying a predicate, not just literal matches. A level specification such as `Infinity` makes `Count` recurse into subexpressions, so it counts every matching leaf at any depth — here `Count[Range[100], _?PrimeQ]` recovers the prime-counting value 25, and the last example tallies how many decimal digits of `2^100` exceed 5.
