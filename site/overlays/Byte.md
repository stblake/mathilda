### Worked examples

```mathematica
In[1]:= Byte  (* the type token stays inert — it has no value of its own *)
```

```mathematica
In[1]:= Head[Byte]  (* it is a plain Symbol *)
```

```mathematica
In[1]:= Attributes[Byte]  (* and carries no attributes *)
```

### Notes

`Byte` is a `Read`/`ReadList` type specification: it reads a single byte, returned as an integer code from 0 to 255. A read that uses it, such as
`ReadList["data.txt", Byte]`, needs a file (or open stream) on disk, so the runnable examples above show
only that the token itself is an inert symbol — it returns raw byte codes (so the character "A" reads as 65), unlike `Number`, which parses numeric tokens.
