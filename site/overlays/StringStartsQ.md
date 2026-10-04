### Worked examples

```mathematica
In[1]:= StringStartsQ["https://x", "http"]  (* a matching prefix *)
```

```mathematica
In[1]:= StringStartsQ["ftp://x", "http"]  (* no match *)
```

```mathematica
In[1]:= StringStartsQ[{"abc", "xyz"}, "a"]  (* threads over a list *)
```

### Notes

`StringStartsQ[s, patt]` is equivalent to
`StringContainsQ[s, StartOfString ~~ patt]`: it anchors the pattern to the start
of the string by wrapping it as `StringExpression[StartOfString, patt]` before
matching.

It shares the predicate core with `StringContainsQ`, `StringFreeQ`, and
`StringEndsQ`, takes an `IgnoreCase` option, and offers the curried operator form
`StringStartsQ[patt]`.
