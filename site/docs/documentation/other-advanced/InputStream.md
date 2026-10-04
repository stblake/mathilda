# InputStream

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`InputStream["name", n]`**

is an object representing an open input stream, as returned by OpenRead. The integer n is an internal handle into the stream registry.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

OpenRead returns an InputStream

```mathematica
In[1]:= s = OpenWrite["/tmp/mathilda_in.txt"]; WriteString[s, "x"]; Close[s]; in = OpenRead["/tmp/mathilda_in.txt"]; h = Head[in]; Close[in]; h
Out[1]= InputStream
```

"name" is the source, 3 the registry handle

```mathematica
In[2]:= FullForm[InputStream["name", 3]]
Out[2]= InputStream["name", 3]
```

## Implementation notes

**Definition.** `InputStream["name", n]` is the object that represents an open input
stream, as returned by `OpenRead`. It has no builtin and no value of its own — it is an
inert stream handle. `"name"` is the source (a file name), and the integer `n` is an
internal handle into the stream registry.

**Representation.** A two-argument `EXPR_FUNCTION` with head `InputStream` (interned
`SYM_InputStream`). `src/io/streams.c` builds it from an entry in the open-stream
registry (`InputStream[name, id]` for an input stream, `OutputStream[name, id]` for an
output stream) and resolves it back to that registry entry whenever `Read`, `ReadList`,
`SetStreamPosition`, `Skip` or `Close` is handed a stream. `Read`/`ReadList` also accept
a bare `"file"` or `File["file"]`, which they open, read, and leave open as an
`InputStream`.

**Usage & limits.** The handle is valid only while the underlying stream is open;
`Close` invalidates it. Reading from an already-open `InputStream` continues from its
current point and leaves it open (unlike a named-file read, which `ReadList` opens and
closes itself). The companion output object is `OutputStream`.

**Attributes:** `Protected`.

## References

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`InputStream["name", n]` is the object representing an open input stream, as returned by
`OpenRead`. `"name"` is the source (a file name) and the integer `n` is an internal
handle into the stream registry. It is an inert handle — no value, no evaluation rule.

The handle is valid only while the stream is open; `Close` invalidates it. Reading from
an already-open `InputStream` continues from its current point and leaves it open, unlike
a named-file read (`Read["file", ...]`), which opens and closes the file itself. The
companion output object is `OutputStream`.
