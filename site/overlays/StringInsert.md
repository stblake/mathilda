### Worked examples

```mathematica
In[1]:= StringInsert["abcdef", "-", 3]  (* insert before the 3rd character *)
```

```mathematica
In[1]:= StringInsert["abcdef", "-", {2, 4}]  (* a copy at each position *)
```

```mathematica
In[1]:= StringInsert["abc", "X", -1]  (* a negative position counts from the end *)
```

### Notes

`StringInsert[s, new, n]` makes the first character of `new` the n-th character
of the result (i.e. it inserts *before* original character n), so a positive `n`
maps to offset `n - 1` and a negative `-k` to `len + n + 1`.

All positions refer to the *original* string, before any insertion; they are
resolved up front into a per-position count, so every copy lands relative to the
untouched input. Valid positions run `1` to `len + 1` (and the negative mirror).
