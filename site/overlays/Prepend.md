### Worked examples

```mathematica
In[1]:= Prepend[{2, 3, 4}, 1]
In[2]:= Prepend[f[b, c], a]  (* works on any head, not just List *)
In[3]:= Prepend[<|a -> 1, b -> 2|>, c -> 3]  (* a rule goes to the front of the association *)
```

### Notes

`Prepend` adds the element at the front while keeping the original head, so it is not
restricted to lists. For an association the argument must be a rule (or rules); it is
placed first, and if its key already exists the old entry is dropped from its former
position. `Append` is the tail-side counterpart.
