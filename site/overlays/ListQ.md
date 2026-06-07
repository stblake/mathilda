### Worked examples

```mathematica
In[1]:= ListQ[{1, 2, 3}]
Out[1]= True

In[2]:= ListQ[5]
Out[2]= False
```

### Notes

`ListQ` tests only whether the head is `List`; it does not inspect the elements.
