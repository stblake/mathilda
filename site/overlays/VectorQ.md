### Worked examples

```mathematica
In[1]:= VectorQ[{1, 2, 3}]
Out[1]= True

In[2]:= VectorQ[{{1, 2}, {3, 4}}]
Out[2]= False
```

### Notes

`VectorQ` is `True` for a flat list none of whose elements are themselves lists; a list of lists (a matrix) is not a vector.
