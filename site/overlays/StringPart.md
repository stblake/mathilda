### Worked examples

```mathematica
In[1]:= StringPart["hello", 1]  (* the first character *)
```

```mathematica
In[1]:= StringPart["hello", -1]  (* the last character *)
```

```mathematica
In[1]:= StringPart["hello", {1, 3, 5}]  (* a list of positions *)
```

### Notes

`StringPart` indexes a string by byte, 1-based, with negative indices counting
from the end (`len + k + 1`). A single index gives a length-1 string; a list of
indices gives a list of characters, and a `Span[m, n, s]` gives the stepped
characters.

A first argument that is itself a list of strings is handled per element. An
out-of-range or non-integer index leaves the call unevaluated.
