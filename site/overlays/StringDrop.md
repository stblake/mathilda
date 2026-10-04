### Worked examples

```mathematica
In[1]:= StringDrop["abcdef", 2]  (* drop the first two characters *)
```

```mathematica
In[1]:= StringDrop["abcdef", -2]  (* drop the last two *)
```

```mathematica
In[1]:= StringDrop["abcdef", {2, 4}]  (* drop a range *)
```

### Notes

`StringDrop` is the complement of `StringTake`: the same sequence specification
selects the characters to *remove*, and the survivors are concatenated in order.
Internally it marks a keep-mask over the bytes and rebuilds the kept ones.

Negative indices count from the end; indexing is byte-based. A decreasing range
`{m, n}` with `m > n` is empty and drops nothing; an out-of-range position leaves
the call unevaluated.
