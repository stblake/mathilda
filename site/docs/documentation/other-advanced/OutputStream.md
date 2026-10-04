# OutputStream

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`OutputStream["name", n]`**

is an object representing an open output stream, as returned by OpenWrite or OpenAppend. The integer n is an internal handle into the stream registry.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

A stream object stays inert — it is a handle, not data

```mathematica
In[1]:= OutputStream["stdout", 1]
Out[1]= OutputStream["stdout", 1]
```

Its head is OutputStream

```mathematica
In[2]:= Head[OutputStream["file.txt", 3]]
Out[2]= OutputStream
```

A protected stream-object head

```mathematica
In[3]:= MemberQ[Attributes[OutputStream], Protected]
Out[3]= True
```

## Implementation notes

**Definition.** `OutputStream["name", n]` is the inert object representing an open
output stream, as returned by `OpenWrite` or `OpenAppend`. The string is the
stream's name (a file name, or `"stdout"`/`"stderr"`); the integer `n` is an
internal handle into the process-global stream registry. It is `Protected`, has no
builtin of its own, and its docstring lives in `info.c`.

**Representation.** `src/io/streams.c` keeps a process-global registry of open
streams. Each open stream is handed to the language as an `OutputStream[name, id]`
object (its input counterpart is `InputStream[name, id]`), where `id` addresses the
registry slot holding the live `FILE*`. The write and stream builtins (`Write`,
`WriteString`, `Close`, ...) resolve a target by looking that `id` up in the
registry, matching the head against the interned `SYM_OutputStream`. As a
two-argument expression its `Head` is `OutputStream`.

**Usage & limits.** It is a handle, not data: you pass it to the write and stream
builtins, and `Close` frees its registry slot. The integer handle is an internal
detail, meaningful only to the registry within the same session. Constructing one by
hand does **not** open a stream — a live object comes from `OpenWrite`/`OpenAppend`.

**Attributes:** `Protected`.

## References

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`OutputStream["name", n]` is the object representing an open output stream: the
string names it (a file, or `"stdout"`/`"stderr"`) and the integer `n` is an
internal handle into the process-global stream registry. A *live* object is produced
by `OpenWrite` or `OpenAppend` (which need a writable file), then passed to `Write`,
`WriteString`, or `Close`; the examples above therefore show only the inert object
itself. The integer handle is an internal registry detail — constructing an
`OutputStream[...]` by hand does not open anything.
