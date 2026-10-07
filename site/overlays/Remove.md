### Worked examples

```mathematica
In[1]:= tmpvar = 7
In[2]:= Remove[tmpvar]
In[3]:= tmpvar
```

### Notes

`Remove[s]` deletes the symbol `s` from the symbol table entirely — a stronger
erase than `ClearAll`, which empties a symbol but keeps its entry. After removal
the name no longer exists; the next reference to it (as in `In[3]`) creates a
fresh, undefined symbol. The result is `Null`.

Arguments may be symbols, strings, or a flat list of them. A `Protected`
symbol is skipped, so `Remove` can never delete a built-in; `Remove`
itself is `Protected`.
