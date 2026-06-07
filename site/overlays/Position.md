### Worked examples

```mathematica
In[1]:= Position[{a,b,a},a]
Out[1]= {{1}, {3}}

In[2]:= Position[{1,2,3,4},_?EvenQ]
Out[2]= {{2}, {4}}
```

### Notes

Returns a list of position specifications (each itself a list), suitable for use with `Extract`.
