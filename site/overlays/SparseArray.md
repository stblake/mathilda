### Worked examples

```mathematica
In[1]:= SparseArray[{1 -> a, 3 -> c, 5 -> e}, 5]  (* stays inert: there is no sparse storage *)
In[2]:= Head[SparseArray[{1 -> 5}, 3]]  (* confirming it does not evaluate *)
In[3]:= Normal[SparseArray[{1 -> a, 3 -> c, 5 -> e}, 5]]  (* Normal expands to the dense list *)
In[4]:= Normal[SparseArray[{{1, 1} -> 1, {2, 2} -> 1}, {2, 2}]]  (* a 2x2 identity *)
In[5]:= Normal[SparseArray[{i_, i_} -> 1, {3, 3}]]  (* pattern rule: the 3x3 identity *)
In[6]:= Normal[SparseArray[{2 -> 7}, 4, -1]]  (* a non-zero default *)
In[7]:= Normal[SparseArray[Band[{1, 1}] -> 1, {3, 3}]]  (* Band fills a diagonal *)
```

### Notes

Mathilda has **no sparse storage**: `SparseArray[…]` is an inert symbolic
specification that stays unevaluated at the REPL. To get values out, wrap it in
`Normal`, which builds the dense nested `List` the specification denotes. Because the
object never becomes a compressed buffer, arithmetic and linear algebra are not
accelerated on it — convert with `Normal` first. `Normal` declines (and leaves the
call unevaluated) when the dense array would exceed 2^27 elements or rank 32.
