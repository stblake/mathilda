### Worked examples

```mathematica
In[1]:= Count[{1,2,1,3,1},1]
Out[1]= 3

In[2]:= Count[{1,2,3,4,5,6},_?EvenQ]
Out[2]= 3

In[3]:= Count[{a,b,a,c,a},a]
Out[3]= 3
```

### Notes

The second argument is a pattern, so `Count[list, _?EvenQ]` counts elements satisfying a predicate, not just literal matches.
