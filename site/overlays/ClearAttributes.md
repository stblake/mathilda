### Worked examples

```mathematica
In[1]:= SetAttributes[g, Orderless]  (* give g an attribute to remove *)
In[1]:= Attributes[g]
In[1]:= ClearAttributes[g, Orderless]  (* now take it away *)
In[1]:= Attributes[g]
```

### Notes

`ClearAttributes[sym, attr]` removes the named attribute bitflags from `sym`; the second
argument may be a single attribute or a list of them, and the first may be one symbol or a
list of symbols. It returns `Null`, so the two `Attributes[g]` queries above show the
before and after.

`ClearAttributes` holds its first argument (`HoldFirst`), so the symbol is cleared rather
than its value. It is the inverse of `SetAttributes`.
