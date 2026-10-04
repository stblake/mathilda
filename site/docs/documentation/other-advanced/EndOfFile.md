# EndOfFile

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`EndOfFile`**

is the symbol returned by Read at the end of a file. ReadList uses it to fill the unread slots of a type sequence truncated by end of file.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

A bare sentinel symbol, not a function

```mathematica
In[1]:= Head[EndOfFile]
Out[1]= Symbol
```

```mathematica
In[2]:= str = OpenWrite["/tmp/mathilda_eof.txt"]; WriteString[str, "only\n"]; Close[str];

In[3]:= ins = OpenRead["/tmp/mathilda_eof.txt"]; Read[ins, Word]; eof = Read[ins, Word]; Close[ins]; eof
Out[3]= EndOfFile
```

A plain marker value, compared by identity

```mathematica
In[4]:= EndOfFile === EndOfFile
Out[4]= True
```

## Implementation notes

**Definition.** `EndOfFile` is the **sentinel symbol** that `Read` returns when a stream
has no more to read. It is a bare symbol, not a function — no builtin, no DownValues, not
even the `Protected` attribute — carrying only the docstring set in `info_init`
(`src/info.c`). It is produced, not consumed, by the reader `src/io/read.c`.

**Representation.** `EndOfFile` stays an inert `EXPR_SYMBOL`; evaluated on its own it is
just `EndOfFile`, and `Head[EndOfFile]` is `Symbol`. The reader constructs it with
`expr_new_symbol(SYM_EndOfFile)` at the two points where a read runs out of input (the
scalar case, and the fallthrough that ends a `ReadList` loop).

**Usage & limits.** When `Read[stream, type]` is called at or past end of file it returns
`EndOfFile`; `ReadList` uses the same sentinel to fill the unread trailing slots of a type
sequence that end-of-file truncated, and stops its loop when it is produced. Code reading a
stream incrementally tests for it to know when to stop — `While[(x = Read[s, Word]) =!=
EndOfFile, ...]`. It is a plain marker value with no numeric or structural meaning of its
own.

**Attributes:** none registered.

## References

- Source: [`src/io/read.c`](https://github.com/stblake/mathilda/blob/main/src/io/read.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`EndOfFile` is the sentinel symbol `Read` returns at the end of a stream. It is an inert
marker — no builtin, no DownValues — *produced* by the reader (`src/io/read.c`), not
consumed: the second `Read[ins, Word]` above has already consumed the only word, so it
hands back `EndOfFile`.

`ReadList` uses the same sentinel to pad the trailing slots of a type sequence that
end-of-file truncated, and stops its loop when it appears. Incremental readers test for it
by identity, as in `While[(x = Read[s, Word]) =!= EndOfFile, ...]`. It carries no numeric or
structural meaning of its own.
