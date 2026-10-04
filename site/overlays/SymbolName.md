### Worked examples

```mathematica
In[1]:= SymbolName[xyz]  (* the symbol's own name as a string *)
```

```mathematica
In[1]:= SymbolName[Global`foo]  (* the context prefix is stripped *)
```

```mathematica
In[1]:= Head[SymbolName[abc]]  (* the result is a String *)
```

### Notes

`SymbolName[sym]` gives the short name of `sym` as a string, with any context
prefix (`Global\``, `P\`Private\``, ...) removed — only the part after the last
backtick is kept. The argument is evaluated first, so to name a symbol that has a
value, pass it held or use the bare symbol. A non-symbol argument is left
unevaluated.
