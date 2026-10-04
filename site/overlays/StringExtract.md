### Worked examples

```mathematica
In[1]:= StringExtract["alpha beta gamma", 2]  (* the 2nd whitespace-delimited block *)
```

```mathematica
In[1]:= StringExtract["a-b-c-d", "-" -> 3]  (* split on "-", take the 3rd block *)
```

```mathematica
In[1]:= StringExtract["one two three four", 2 ;; 3]  (* a span of blocks *)
```

### Notes

`StringExtract` splits a string into blocks and selects by position. It is built
directly on `StringSplit` (it synthesises and evaluates `StringSplit[str, sep]`),
so the split behaviour is identical and `StringExtract[s, sep -> All]` is exactly
`StringSplit[s, sep]`.

A bare position gets a depth-default separator (whitespace at the lowest level,
growing runs of `"\n"` above). An out-of-range index yields
`Missing["PartAbsent", n]` rather than leaving the call unevaluated.
