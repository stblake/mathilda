### Worked examples

```mathematica
In[1]:= StringRepeat["ab", 4]  (* four copies *)
```

```mathematica
In[1]:= StringRepeat["-", 10]  (* a rule of ten dashes *)
```

```mathematica
In[1]:= StringRepeat["abc", 5, 7]  (* capped at 7 bytes, so the last copy is partial *)
```

### Notes

`StringRepeat[s, n]` concatenates `n` copies of `s`. The optional third argument
caps the total length at `max` bytes, keeping a partial final copy; the buffer is
filled cyclically (`buf[i] = s[i mod len]`) so truncation falls out naturally.

`n == 0`, an empty base string, or `max == 0` gives `""`. The `n * len` product is
guarded against overflow, which is tolerated only when a `max` cap keeps the
result finite.
