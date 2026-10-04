### Worked examples

```mathematica
In[1]:= Number  (* the type token stays inert — it has no value of its own *)
```

```mathematica
In[1]:= Head[Number]  (* it is a plain Symbol *)
```

```mathematica
In[1]:= Attributes[Number]  (* and carries no attributes *)
```

### Notes

`Number` is a `Read`/`ReadList` type specification: it reads a number, returned as an integer when the token has no decimal point or exponent and as an approximate number otherwise. A read that uses it, such as
`ReadList["data.txt", Number]`, needs a file (or open stream) on disk, so the runnable examples above show
only that the token itself is an inert symbol — an integer-looking token reads back as an `Integer` and a decimal or exponent token as a `Real`.
