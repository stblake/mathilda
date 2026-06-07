### Worked examples

```mathematica
In[1]:= $Epilog := Print["bye"]
Out[1]= Null

In[2]:= Quit[]
"bye"
```

### Notes

`$Epilog` is a symbol whose value, if any, is evaluated exactly once when the
session terminates via `Quit[]` or EOF. Assigning it with `:=` (delayed) lets
it run cleanup or a parting action at shutdown; the `"bye"` above is printed as
the session exits.
