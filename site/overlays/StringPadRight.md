### Worked examples

```mathematica
In[1]:= StringPadRight["42", 5, "0"]  (* right-pad to width 5 with zeros *)
```

```mathematica
In[1]:= StringPadRight["abc", 7, ".-"]  (* a multi-character pad cycles *)
```

```mathematica
In[1]:= StringPadRight["abcdef", 3]  (* too long: truncation keeps the FIRST n *)
```

### Notes

`StringPadRight` is the mirror of `StringPadLeft`: it pads on the right and, when
truncating, keeps the first `n` bytes. The two share one implementation
(`pad_dispatch`) selected by a `left` flag.

The pad string is laid down cyclically from the left, so a multi-character pad
such as `".-"` repeats `.-.-...`. Lengths are byte counts.
