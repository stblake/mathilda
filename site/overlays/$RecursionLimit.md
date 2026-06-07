### Worked examples

```mathematica
In[1]:= $RecursionLimit
Out[1]= 1024

In[2]:= $RecursionLimit = 500
Out[2]= 500
```

### Notes

`$RecursionLimit` gives the maximum depth of nested evaluator invocations;
its default is `1024`. Assigning a positive integer of at least `20` updates
the limit, while smaller values are rejected with a `$RecursionLimit::limset`
message.
