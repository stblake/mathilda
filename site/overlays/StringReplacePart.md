### Worked examples

```mathematica
In[1]:= StringReplacePart["abcdef", "XY", {2, 3}]  (* replace a character range *)
```

```mathematica
In[1]:= StringReplacePart["abcdef", "", {2, 4}]  (* an empty string deletes the range *)
```

```mathematica
In[1]:= StringReplacePart["abcdefgh", {"1", "2"}, {{1, 2}, {5, 6}}]  (* one new string per range *)
```

### Notes

Position specifications are `{m, n}` first/last character pairs — the form
`StringPosition` returns — with negative positions counting from the end. All
positions refer to the *original* string, before any replacement.

A single new string is broadcast to every range; a list of new strings must match
the number of ranges. Overlapping ranges are not allowed: a later range touching
an already-claimed position triggers `StringReplacePart::ovlp` and is dropped. The
operator form is `StringReplacePart[new, part][old]`.
