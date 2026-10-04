# Word

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Word`**

is a type specification in Read and ReadList that reads a sequence of characters delimited by word separators (see WordSeparators).

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

A bare read-type token, not a function

```mathematica
In[1]:= Head[Word]
Out[1]= Symbol
```

```mathematica
In[2]:= str = OpenWrite["/tmp/mathilda_word.txt"]; WriteString[str, "alpha beta gamma\n"]; Close[str];
```

Each whitespace-delimited run as a string

```mathematica
In[3]:= ReadList["/tmp/mathilda_word.txt", Word]
Out[3]= {"alpha", "beta", "gamma"}
```

Just the first word

```mathematica
In[4]:= ins = OpenRead["/tmp/mathilda_word.txt"]; w = Read[ins, Word]; Close[ins]; w
Out[4]= "alpha"
```

## Implementation notes

**Definition.** `Word` is a **read type** token for `Read` and `ReadList`: it asks the
reader to return the next run of characters delimited by word separators as a string. It is
a bare symbol, not a function — no builtin, no DownValues, not even the `Protected`
attribute — carrying only the docstring set in `info_init` (`src/info.c`). Its meaning
lives in the reader `src/io/read.c`, where `read_type_of` maps the symbol `SYM_Word` to the
internal `RT_WORD` case.

**Representation.** `Word` stays an inert `EXPR_SYMBOL`; evaluated on its own it is just
`Word`, and `Head[Word]` is `Symbol`. It appears only as a type argument, e.g.
`Read[stream, Word]` or `ReadList["file", Word]`, and inside the type *structure* of a
nested read such as `ReadList["file", {Word, Number}]`.

**Usage & limits.** When the reader hits a `Word` slot it skips leading word separators
(by default space and tab, configurable through `WordSeparators`), then collects characters
up to the next separator and returns them as a `String`. `TokenWords` can force specific
strings to read as standalone words, and `NullWords -> True` keeps empty fields between
adjacent separators. A `Word` read past end of file returns `EndOfFile`. `Word` is distinct
from `Record` (a whole record up to a record separator) and `String` (a line).

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Word` is a **read type** for `Read` and `ReadList`: it reads the next run of characters
delimited by word separators and returns it as a string. It is an inert token symbol — no
builtin, no DownValues — recognised by the reader in `src/io/read.c`, which maps it to its
internal word case.

`ReadList["file", Word]` returns every word as a flat list of strings; `Read[stream, Word]`
returns just the next one and advances the stream. The separators default to space and tab
and are set by `WordSeparators`; `TokenWords` forces given strings to read as standalone
words, and `NullWords -> True` keeps empty fields between adjacent separators. A `Word` read
past end of file returns `EndOfFile`.
