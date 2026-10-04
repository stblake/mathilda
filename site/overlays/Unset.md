### Worked examples

```mathematica
In[1]:= v = 10
In[2]:= v =.
In[3]:= v
```

```mathematica
In[1]:= g[1] = a; g[2] = b; g[1] =.; DownValues[g]  (* removes only the g[1] rule *)
```

### Notes

`lhs =.` (`Unset`) removes the one definition whose left-hand side is `lhs`,
rather than every rule on a symbol. For a bare symbol `v =.` drops its OwnValue;
for a pattern `g[1] =.` drops exactly that DownValue, leaving the others in place
— here `DownValues[g]` keeps only the `g[2]` rule. The result is `Null`.

An unassignable left-hand side is left alone, and a `Protected` or `Locked`
symbol is refused with `Unset::wrsym`. `Unset` is `HoldFirst`, so the target is
not evaluated to its value before the matching rule is located.
