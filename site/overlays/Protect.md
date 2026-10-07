### Worked examples

```mathematica
In[1]:= Protect[foo, bar]  (* returns the names whose Protected attribute was newly set *)
In[2]:= MemberQ[Attributes[foo], Protected]
```

### Notes

`Protect[s1, s2, ...]` sets the `Protected` attribute on each symbol, which is
what the evaluator and `Set` consult to refuse redefinition. It returns a list of
the names (as strings) whose state actually *changed*, so re-protecting an
already-protected symbol gives `{}`.

Arguments may be symbols, strings, or a flat list of them. `Protected` is the
attribute every built-in carries; `Unprotect`
is its inverse and the usual first step before extending a built-in's behaviour.
