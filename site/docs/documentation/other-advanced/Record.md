# Record

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Record`**

is a type specification in Read and ReadList that reads a sequence of characters delimited by record separators (see RecordSeparators).

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

The type token stays inert — it has no value of its own

```mathematica
In[1]:= Record
Out[1]= Record
```

It is a plain Symbol

```mathematica
In[2]:= Head[Record]
Out[2]= Symbol
```

And carries no attributes

```mathematica
In[3]:= Attributes[Record]
Out[3]= {}
```

## Implementation notes

**Definition.** `Record` is a type specification used in `Read` and `ReadList`: it reads a sequence of characters delimited by record separators (see `RecordSeparators`), returned as a string.
It is an inert symbol — no builtin, no value, no attributes — declared with its
docstring in `info.c` and recognised by the reader.

**Representation.** A bare `EXPR_SYMBOL`. The read-type scanner in `src/io/read.c`
maps it to an internal read-type tag (`if (n == SYM_Record) *out = RT_...`), after
which the matching reader pulls the value from the stream. The record boundaries come from the `RecordSeparators` option, not from `Byte`/`Number` parsing. As a plain symbol
its `Head` is `Symbol` and it carries no attributes.

**Usage & limits.** Meaningful only as the type argument of `Read`/`ReadList` over a
file or an open stream; standing alone it is just a symbol. each record is returned as a `String`, split on the `RecordSeparators` strings.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Record` is a `Read`/`ReadList` type specification: it reads a sequence of characters delimited by record separators (see `RecordSeparators`), returned as a string. A read that uses it, such as
`ReadList["data.txt", Record, RecordSeparators -> {","}]`, needs a file (or open stream) on disk, so the runnable examples above show
only that the token itself is an inert symbol — each record is returned as a `String`, split on the `RecordSeparators` strings.
