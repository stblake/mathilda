### Worked examples

```mathematica
In[1]:= Flatten[{{1,2},{3,{4}}}]
Out[1]= {1, 2, 3, 4}

In[2]:= Flatten[{{1,{2,3}},{4,{5}}},1]
Out[2]= {1, {2, 3}, 4, {5}}
```

### Notes

Without a level argument, `Flatten` removes all nesting; with a level `n` it flattens only the top `n` levels.
