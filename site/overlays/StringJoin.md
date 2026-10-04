### Worked examples

```mathematica
In[1]:= StringJoin["abc", "def"]  (* concatenate the arguments *)
```

```mathematica
In[1]:= "foo" <> "bar" <> "baz"  (* the infix <> operator *)
```

```mathematica
In[1]:= StringJoin[{"a", "b", "c"}]  (* nested lists are flattened *)
```

### Notes

`StringJoin` gathers every leaf string (descending through any `List` wrappers),
sums their lengths, and copies them into one buffer. It is `Flat` and
`OneIdentity`, so the evaluator flattens nested `StringJoin` and `<>` before the
builtin runs.

The zero-argument form `StringJoin[]` is `""`. Any non-string, non-`List` leaf
leaves the call unevaluated.
