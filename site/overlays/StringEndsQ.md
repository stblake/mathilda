### Worked examples

```mathematica
In[1]:= StringEndsQ["report.pdf", ".pdf"]  (* a matching suffix *)
```

```mathematica
In[1]:= StringEndsQ["report.txt", ".pdf"]  (* no match *)
```

```mathematica
In[1]:= StringEndsQ[{"a.c", "b.h", "d.c"}, ".c"]  (* threads over a list *)
```

### Notes

`StringEndsQ[s, patt]` is equivalent to
`StringContainsQ[s, patt ~~ EndOfString]`: it anchors the pattern to the end of
the string by wrapping it as `StringExpression[patt, EndOfString]` before
matching.

It shares the predicate core with `StringContainsQ`, `StringFreeQ`, and
`StringStartsQ`, takes an `IgnoreCase` option, and offers the curried operator
form `StringEndsQ[patt]`.
