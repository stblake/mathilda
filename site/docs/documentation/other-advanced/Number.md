# Number

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Number`**

is a type specification in Read and ReadList that reads a number, returned as an integer when the token has no decimal point or exponent and as an approximate number otherwise.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

The type token stays inert — it has no value of its own

```mathematica
In[1]:= Number
Out[1]= Number
```

It is a plain Symbol

```mathematica
In[2]:= Head[Number]
Out[2]= Symbol
```

And carries no attributes

```mathematica
In[3]:= Attributes[Number]
Out[3]= {}
```

## Implementation notes

**Definition.** `Number` is a type specification used in `Read` and `ReadList`: it reads a number, returned as an integer when the token has no decimal point or exponent and as an approximate number otherwise.
It is an inert symbol — no builtin, no value, no attributes — declared with its
docstring in `info.c` and recognised by the reader.

**Representation.** A bare `EXPR_SYMBOL`. The read-type scanner in `src/io/read.c`
maps it to an internal read-type tag (`if (n == SYM_Number) *out = RT_...`), after
which the matching reader pulls the value from the stream.  As a plain symbol
its `Head` is `Symbol` and it carries no attributes.

**Usage & limits.** Meaningful only as the type argument of `Read`/`ReadList` over a
file or an open stream; standing alone it is just a symbol. an integer-looking token reads back as an `Integer` and a decimal or exponent token as a `Real`.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Number` is a `Read`/`ReadList` type specification: it reads a number, returned as an integer when the token has no decimal point or exponent and as an approximate number otherwise. A read that uses it, such as
`ReadList["data.txt", Number]`, needs a file (or open stream) on disk, so the runnable examples above show
only that the token itself is an inert symbol — an integer-looking token reads back as an `Integer` and a decimal or exponent token as a `Real`.
