### Worked examples

```mathematica
In[1]:= StringFreeQ["abcdef", "xyz"]  (* nothing matches, so free *)
```

```mathematica
In[1]:= StringFreeQ["abcdef", "cd"]  (* a substring matches, so not free *)
```

```mathematica
In[1]:= StringFreeQ["hello", LetterCharacter]  (* letters are present *)
```

### Notes

`StringFreeQ` is the exact negation of `StringContainsQ` and shares its
implementation (one `SqKind` tag selects the inverting answer). It is `True` when
no substring matches the pattern.

It takes an `IgnoreCase` option and offers the curried operator form
`StringFreeQ[patt]`. A list of patterns means "free of any of them".
