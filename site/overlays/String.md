### Worked examples

```mathematica
In[1]:= Head["abc"]  (* String is the head of string objects *)
```

```mathematica
In[1]:= MatchQ["abc", _String]  (* so _String matches any string *)
```

```mathematica
In[1]:= Head[String]  (* the bare symbol String is itself a Symbol *)
```

### Notes

`String` is the head carried by string leaves — `Head["text"]` is `String`, which
is what makes `_String` match any string in a pattern, type test, or
`Cases`/`Select`. Strings are stored as `EXPR_STRING` nodes, not as `String[...]`
compounds, so `String` appears as a head but never wraps anything. It is *also* a
`Read`/`ReadList` type specification that reads one line (up to a newline); that
role needs a file or open stream to exercise. The bare symbol `String` is an
ordinary `Symbol`.
