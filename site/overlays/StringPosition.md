### Worked examples

```mathematica
In[1]:= StringPosition["abcabc", "bc"]  (* {start, end} pairs, 1-based inclusive *)
```

```mathematica
In[1]:= StringPosition["aaaa", "aa"]  (* Overlaps -> True is the default here *)
```

```mathematica
In[1]:= StringPosition["a1b2c3", DigitCharacter]  (* a character-class pattern *)
```

### Notes

`StringPosition` returns `{start, end}` character positions in the 1-based
inclusive form that `StringTake`, `StringDrop`, and `StringReplacePart` consume.
It shares the `regex_scan` enumerator with `StringCases`/`StringCount`.

Unlike those two it defaults to `Overlaps -> True` (matching the Wolfram
Language), so overlapping matches at distinct starts are listed. An optional
third integer argument keeps only the first `n` matches.
