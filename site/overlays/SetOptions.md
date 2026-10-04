### Worked examples

```mathematica
In[1]:= Options[FactorInteger]
In[2]:= SetOptions[FactorInteger, GaussianIntegers -> True]
In[3]:= Options[FactorInteger]
```

### Notes

`SetOptions[s, name -> value, ...]` changes a symbol's *default* option settings
and returns the full updated option list. The change is made in place: the
matching option keeps its original position and only its value is replaced, so
later `Options[s]` and option-reading builtins see the new default.

`SetOptions` can only change options a symbol already has — a name that is not a
known option raises `SetOptions::optnf` and leaves the settings untouched. The
first argument must be a symbol, and a `Locked` symbol is refused with
`SetOptions::locked`.
