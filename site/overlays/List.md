### Worked examples

```mathematica
In[1]:= FullForm[{p, q, r}]  (* {...} is syntax for List[...] *)
```

```mathematica
In[1]:= Head[{1, 2, 3}]
```

```mathematica
In[1]:= List[1, 2, 3] === {1, 2, 3}  (* the long form and the brace form are identical *)
```

### Notes

`List[e1, e2, ...]`, written `{e1, e2, ...}`, is the fundamental ordered container.
Vectors are lists, matrices are lists of lists, and the structural operators (`Part`,
`Map`, `Take`, `Drop`, `Length`, ...) all act on `List`. The parser rewrites `{...}` to
`List[...]` and the printer renders it back, so `FullForm` is where the `List` head
shows through.

Elements are evaluated normally and kept in the given order — `List` has **no**
`Orderless` attribute, so `{3, 1, 2}` stays `{3, 1, 2}`. In fact it carries no
attributes at all (`Attributes[List]` is `{}`). This uniformity is the design: because
everything is an expression, one set of generic tools works on every list.
