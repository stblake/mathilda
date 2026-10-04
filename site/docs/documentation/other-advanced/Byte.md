# Byte

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Byte`**

is a type specification in Read and ReadList that reads a single byte, returned as an integer code from 0 to 255.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

The type token stays inert — it has no value of its own

```mathematica
In[1]:= Byte
Out[1]= Byte
```

It is a plain Symbol

```mathematica
In[2]:= Head[Byte]
Out[2]= Symbol
```

And carries no attributes

```mathematica
In[3]:= Attributes[Byte]
Out[3]= {}
```

## Implementation notes

**Definition.** `Byte` is a type specification used in `Read` and `ReadList`: it reads a single byte, returned as an integer code from 0 to 255.
It is an inert symbol — no builtin, no value, no attributes — declared with its
docstring in `info.c` and recognised by the reader.

**Representation.** A bare `EXPR_SYMBOL`. The read-type scanner in `src/io/read.c`
maps it to an internal read-type tag (`if (n == SYM_Byte) *out = RT_...`), after
which the matching reader pulls the value from the stream. Because `Byte` reads the raw stream, it returns the numeric code of each octet, not a parsed number. As a plain symbol
its `Head` is `Symbol` and it carries no attributes.

**Usage & limits.** Meaningful only as the type argument of `Read`/`ReadList` over a
file or an open stream; standing alone it is just a symbol. it returns raw byte codes (so the character "A" reads as 65), unlike `Number`, which parses numeric tokens.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Byte` is a `Read`/`ReadList` type specification: it reads a single byte, returned as an integer code from 0 to 255. A read that uses it, such as
`ReadList["data.txt", Byte]`, needs a file (or open stream) on disk, so the runnable examples above show
only that the token itself is an inert symbol — it returns raw byte codes (so the character "A" reads as 65), unlike `Number`, which parses numeric tokens.
