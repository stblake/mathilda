# Character

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Character`**

is a type specification in Read and ReadList that reads a single character, returned as a one-character string.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

Reads one character as a string

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_char.txt"]; WriteString[s, "Hi"]; Close[s]; Read["/tmp/mathilda_char.txt", Character]
Out[1]= "H"
```

One of the eight read types

```mathematica
In[2]:= MemberQ[{Byte, Character, Expression, Number, Real, Record, String, Word}, Character]
Out[2]= True
```

## Implementation notes

**Definition.** `Character` is a **type specification** used by `Read` and `ReadList`.
It has no builtin and no value — it is an inert token that tells the reader to consume
exactly one character from the stream and return it as a one-character string.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_Character`). `read.c`'s
`type_from_symbol` maps it to the internal `RT_CHARACTER` reader, one of the eight read
types (`Byte`, `Character`, `Expression`, `Number`, `Real`, `Record`, `String`,
`Word`). It may appear on its own or nested inside a type-structure argument
(`{Character, Character}`, `Hold[...]`, any head), which the reader fills depth-first.

**Usage & limits.** Meaningful only as a `Read`/`ReadList` type argument; evaluated on
its own it simply returns itself. Past end of file the reader returns `EndOfFile`.
Contrast `Byte`, which returns the raw integer code 0-255 rather than a string.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Character` is a type specification for `Read` and `ReadList`: it consumes exactly one
character from the stream and returns it as a one-character string. It has no value of
its own, so on its own it just evaluates to itself.

It may stand alone or nest inside a type structure (`{Character, Character}`,
`Hold[...]`, any head), which the reader fills depth-first. Contrast `Byte`, which
returns the raw integer code 0-255 instead of a string. Past end of file the reader
returns `EndOfFile`.
