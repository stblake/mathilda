### Worked examples

```mathematica
In[1]:= {a, Nothing, b, Nothing, c}  (* every Nothing simply disappears *)
```

```mathematica
In[1]:= {1, 2, Nothing[x, y], 3}  (* a Nothing[...] form is removed too *)
```

```mathematica
In[1]:= Table[If[PrimeQ[i], i, Nothing], {i, 10}]  (* keep only the elements that pass a test *)
```

### Notes

`Nothing` is a symbol that is automatically removed from any list in which it
appears. It is the identity element for list construction, which makes
`Table[If[test, val, Nothing], ...]` the standard idiom for building a list of
just the values for which a test holds — no `Select` or `DeleteCases` pass is
needed. Any `Nothing[...]` form is stripped likewise.

The removal is list-specific: for a non-`List` head, `Nothing` is an ordinary
symbol and is left in place.
