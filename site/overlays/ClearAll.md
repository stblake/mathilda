### Worked examples

```mathematica
In[1]:= SetAttributes[g, Orderless]; g = 5
In[2]:= ClearAll[g]
In[3]:= {g, Attributes[g]}
```

```mathematica
In[1]:= ClearAll[a, b]  (* several symbols in one call *)
```

### Notes

`ClearAll[s]` is the thorough erase: it removes `s`'s OwnValues and DownValues
*and* its attributes and usage message, so the symbol is returned to its
pristine undefined state — compare `Clear[s]`, which drops only the
values. Here `g` loses both its assigned value and the `Orderless` attribute, so
`Attributes[g]` is `{}`.

Arguments may be symbols, strings naming symbols, or a flat list of them
(`ClearAll[{a, b}]`). A `Protected` or `Locked` symbol is skipped, which is what
keeps `ClearAll` from ever gutting a built-in. The result is `Null`.
