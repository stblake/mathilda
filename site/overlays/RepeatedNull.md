### Worked examples

```mathematica
In[1]:= MatchQ[{}, {RepeatedNull[1]}]
Out[1]= True

In[2]:= MatchQ[{1, 1}, {RepeatedNull[1]}]
Out[2]= True
```

### Notes

`RepeatedNull[p]` (postfix `p...`) is like `Repeated` but matches *zero* or more occurrences, so an empty sequence also succeeds (Out[1]). Use it when the repeated part is optional.
