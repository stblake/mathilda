### Worked examples

```mathematica
In[1]:= StringCases["the cat sat", LetterCharacter ..]  (* runs of letters, left to right *)
```

```mathematica
In[1]:= StringCases["a1b2c3", RegularExpression["[a-z](\\d)"] -> "$1"]  (* keep only the captured digit *)
```

```mathematica
In[1]:= StringCases["AAAA", "AA", Overlaps -> True]  (* overlapping matches, one per start *)
```

### Notes

Matching is non-overlapping and greedy by default (`Overlaps -> False`). A
`patt -> rhs` rule rewrites each match, expanding `$0`/`$n`; a bare pattern
returns the matched substring itself.

`StringCases` shares its match enumerator (`regex_scan`) with `StringCount` and
`StringPosition`, so `StringCount[s, p]` always equals
`Length[StringCases[s, p]]` under any option setting. A list of subjects threads
element-wise.
