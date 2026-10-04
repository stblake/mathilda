### Worked examples

```mathematica
In[1]:= Options[NumberForm, DigitBlock]  (* a single named option as a rule *)
```

```mathematica
In[1]:= Options[gsym]  (* a symbol with no registered options has none *)
```

### Notes

`Options[s]` returns a symbol's default option settings as a list of rules
`{name -> value, ...}`; a symbol with none gives `{}`. `Options[s, name]` or
`Options[s, {names}]` selects just the requested rules, in the order asked for.

Applied to a compound expression rather than a symbol, `Options[expr]` returns
only the option rules that appear explicitly among `expr`'s arguments. Name
matching ignores any context prefix, and the returned list is a fresh copy, so
nothing stored is aliased.
