### Worked examples

```mathematica
In[1]:= StringPadLeft["42", 5, "0"]  (* left-pad to width 5 with zeros *)
```

```mathematica
In[1]:= StringPadLeft["abcdef", 3]  (* too long: truncation keeps the LAST n *)
```

```mathematica
In[1]:= StringPadLeft[{"a", "bb", "ccc"}]  (* a list pads to the longest element *)
```

### Notes

`StringPadLeft` makes a string a given length, padding on the left or truncating
it. When the string is longer than `n` it keeps the last `n` bytes; the pad
string (a single space by default) is laid down cyclically, `p[i mod plen]`.

The one-argument list form pads every element to the length of the longest, so
all come out equal width. Lengths are byte counts. A list-valued pad string is
not supported and leaves the call unevaluated.
