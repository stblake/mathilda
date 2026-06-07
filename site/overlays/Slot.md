### Worked examples

```mathematica
In[1]:= (#1 + #2 &)[3, 4]
Out[1]= 7

In[2]:= Map[#^2 &, {1, 2, 3, 4}]
Out[2]= {1, 4, 9, 16}

In[3]:= FullForm[#2]
Out[3]= Slot[2]
```

### Notes

`#` (or `#1`) is `Slot[1]`, the first argument of the enclosing pure function (`&`); `#n` refers to the n-th argument. A bare `#` inside `Map` (Out[2]) receives each list element in turn.
