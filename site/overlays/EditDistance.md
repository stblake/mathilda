### Worked examples

```mathematica
In[1]:= EditDistance["Sunday", "Saturday"]  (* three edits turn one word into the other *)
```

```mathematica
In[1]:= EditDistance[{1, 2, 3, 4}, {1, 3, 4}]  (* works on lists, not only strings *)
```

```mathematica
In[1]:= EditDistance["", "abc"]  (* three insertions from the empty string *)
```

### Notes

`EditDistance[a, b]` is the Levenshtein distance: the fewest single-element
insertions, deletions, and substitutions that turn `a` into `b`. Both arguments
must be strings, or both lists — a mixed pair is left unevaluated. Elements are
compared for equality, so the same function serves strings (character by
character) and lists of arbitrary expressions.

Strings are compared byte by byte, so a multi-byte UTF-8 character counts as
several elements. For the positional "how many slots differ" count on
equal-length sequences, use `HammingDistance` instead.
