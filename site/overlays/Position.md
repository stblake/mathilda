### Worked examples

```mathematica
In[1]:= Position[{a,b,a},a]
Out[1]= {{1}, {3}}

In[2]:= Position[{1,2,3,4},_?EvenQ]
Out[2]= {{2}, {4}}
```

```mathematica
In[1]:= Position[{1, 2, 3, 4, 5, 6}, _?PrimeQ]
Out[1]= {{2}, {3}, {5}}
```

```mathematica
In[1]:= Position[x^2 + y^2 + z^2, _Symbol]
Out[1]= {{0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {3, 0}, {3, 1}}
```

```mathematica
In[1]:= Position[Sin[Cos[x] + Tan[x]], x, Infinity]
Out[1]= {{1, 1, 1}, {1, 2, 1}}
```

### Notes

Returns a list of position specifications (each itself a list), suitable for use with `Extract`.
