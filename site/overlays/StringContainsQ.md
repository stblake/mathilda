### Worked examples

```mathematica
In[1]:= StringContainsQ["abcdef", "cd"]  (* some substring matches *)
```

```mathematica
In[1]:= StringContainsQ["Mathilda", "math", IgnoreCase -> True]  (* case folded *)
```

```mathematica
In[1]:= Select[{"cat", "dog", "cow"}, StringContainsQ["o"]]  (* operator form as a predicate *)
```

### Notes

`StringContainsQ` searches *inside* the string, where `StringMatchQ` anchors to
the whole of it: it is equivalent to `StringMatchQ[s, ___ ~~ patt ~~ ___]` and to
`!StringFreeQ[s, patt]`.

The one-argument operator form `StringContainsQ[patt]` curries into a predicate
(a curried `IgnoreCase -> True` is carried through), which makes it a natural
second argument to `Select`. It shares a single core with `StringFreeQ`,
`StringStartsQ`, and `StringEndsQ`.
