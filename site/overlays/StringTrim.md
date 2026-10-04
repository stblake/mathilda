### Worked examples

```mathematica
In[1]:= StringTrim["   hello   "]  (* default: trim whitespace from both ends *)
```

```mathematica
In[1]:= StringTrim["xxhelloxx", "x"]  (* each end is trimmed to a fixed point *)
```

```mathematica
In[1]:= StringTrim["...data...", "."]  (* a literal pattern *)
```

### Notes

`StringTrim` removes substrings matching the pattern from the *start and end*
only (never the middle), stripping each end repeatedly to a fixed point — so a
single-character pattern still removes a whole run. The default pattern is
`Whitespace`.

Anchoring is one-ended: the front is matched start-relative and the back with a
`\z` anchor over a truncated length, rather than the whole-string `\A...\z` wrap
the other regex builtins use. The pattern accepts the full shared string-pattern
vocabulary.
