### Worked examples

```mathematica
In[1]:= Record  (* the type token stays inert — it has no value of its own *)
```

```mathematica
In[1]:= Head[Record]  (* it is a plain Symbol *)
```

```mathematica
In[1]:= Attributes[Record]  (* and carries no attributes *)
```

### Notes

`Record` is a `Read`/`ReadList` type specification: it reads a sequence of characters delimited by record separators (see `RecordSeparators`), returned as a string. A read that uses it, such as
`ReadList["data.txt", Record, RecordSeparators -> {","}]`, needs a file (or open stream) on disk, so the runnable examples above show
only that the token itself is an inert symbol — each record is returned as a `String`, split on the `RecordSeparators` strings.
