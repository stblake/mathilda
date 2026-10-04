### Worked examples

```mathematica
In[1]:= AssociationQ[<|a -> 1, b -> 2|>]
```

```mathematica
In[1]:= AssociationQ[{a -> 1, b -> 2}]  (* a list of rules is not an association *)
```

### Notes

`AssociationQ[expr]` returns `True` only for a well-formed association — one with
the `Association` head whose every entry is a two-argument `Rule` or
`RuleDelayed`. A bare list of rules such as `{a -> 1, b -> 2}` is a `List`, not an
association, so it answers `False`; likewise a malformed `Association[1, 2]` left
unevaluated by the constructor. Like every `*Q` predicate it always returns a
Boolean, never staying unevaluated.
