### Worked examples

```mathematica
In[1]:= Protect[baz]
In[2]:= Unprotect[baz]  (* returns the names whose Protected attribute was actually cleared *)
In[3]:= MemberQ[Attributes[baz], Protected]
```

### Notes

`Unprotect[s1, s2, ...]` is the inverse of `Protect`: it clears the `Protected`
attribute from each symbol and returns the list of names (as strings) whose state
changed, so unprotecting a symbol that was never protected gives `{}`.

Once a symbol is unprotected, `Set`/`SetDelayed` and the clearing heads accept it
again — the standard way to override or extend a built-in's definition.
Arguments may be symbols, strings, or a flat list; a `Locked` symbol cannot be
unprotected.
