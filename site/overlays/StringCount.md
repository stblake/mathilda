### Worked examples

```mathematica
In[1]:= StringCount["mississippi", "ss"]  (* non-overlapping by default *)
```

```mathematica
In[1]:= StringCount["a1b2c3d4", DigitCharacter]  (* a character-class pattern *)
```

```mathematica
In[1]:= StringCount["banana", "a", Overlaps -> True]  (* count overlapping starts *)
```

### Notes

`StringCount` is the counting-only companion of `StringCases`: it runs the same
`regex_scan` enumeration but records only a small span per match rather than
building each substring, so it is cheaper for a pure count.

The result is exactly `Length[StringCases[...]]` for every pattern and option.
`Overlaps`/`IgnoreCase` behave as in `StringCases`, and a list of subjects
threads, giving one count each.
