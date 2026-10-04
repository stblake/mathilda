### Worked examples

```mathematica
In[1]:= StringQ[Information[Sin]]  (* Information hands back the docstring as a string *)
```

```mathematica
In[1]:= Head[Information[Plus]]  (* its result is a String *)
```

### Notes

`Information[sym]` returns the symbol's docstring — the same text the interactive `?sym`
shortcut prints — as a string, or a `No information available` string when the symbol has
none. Every builtin registers its docstring via `symtab_set_docstring`, which is the store
both `Information` and `?sym` read.

Because the result is the usage text itself (long, multi-line, and specific to the symbol
asked about), the examples above probe it structurally; in the REPL you would simply type
`?Sin`.
