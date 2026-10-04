# Close

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Close[stream]`**

closes an open input or output stream and returns its file name.

**`Close["file"] closes a stream opened for that file. Returns $Failed if`**

the stream is not open.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_io_close.txt"]
Out[1]= OutputStream["/tmp/mathilda_io_close.txt", 1]
```

The stream is open

```mathematica
In[2]:= Streams["/tmp/mathilda_io_close.txt"]
Out[2]= {OutputStream["/tmp/mathilda_io_close.txt", 1]}
```

Close returns the file name

```mathematica
In[3]:= Close[str]
Out[3]= "/tmp/mathilda_io_close.txt"
```

And removes it from the registry

```mathematica
In[4]:= Streams["/tmp/mathilda_io_close.txt"]
Out[4]= {}
```

### Applications (4)

Open a stream

```mathematica
In[5]:= str = OpenWrite["/tmp/mathilda_close.txt"];
```

Close returns the file name, not Null

```mathematica
In[6]:= Close[str]
Out[6]= "/tmp/mathilda_close.txt"
```

Open, discarding the object

```mathematica
In[7]:= OpenWrite["/tmp/mathilda_close2.txt"];
```

Close by file name instead of by object

```mathematica
In[8]:= Close["/tmp/mathilda_close2.txt"]
Out[8]= "/tmp/mathilda_close2.txt"
```

## Implementation notes

**Algorithm.** `builtin_close` resolves its argument two ways: an `InputStream`/`OutputStream`
object through `stream_handle` (reading the id out of the object's second slot), or a `"file"`
/`File["file"]` name through `stream_find_by_name(name, -1)`, which returns the
**most-recently-opened** matching stream. With the slot found, it captures the file name into a
fresh string, calls `stream_close_id` — which `free`s the slot's `name` and `buf`, `fclose`s
the `fp` if any, and `memset`s the slot so its `id` returns to `0` (free) — and **returns the
file name**. A stream that is not open prints `Close::stream` and returns `$Failed`.

**Data structures — the stream layer.** Operates on the process-global `Stream` registry
described under `OpenRead`. Closing is what frees a slot deterministically; the
monotonic `g_next_id` guarantees the freed id is never handed to a later open, so a stale
object cannot alias a new stream. Returning the name (rather than `Null`) is what lets `Close`
be the convenient last line of a round trip.

**Complexity / limits.** `O(n)` over the registry to find a slot by id or name (`n` = open
streams, tiny in practice). For an output stream the `fclose` is what guarantees the file is
flushed to disk. Pure ANSI C99 (`fclose`). `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## References

**See also:** [InputStream](../../other-advanced/InputStream/), [OutputStream](../../other-advanced/OutputStream/)

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`Close[stream]` closes an open `InputStream`/`OutputStream` and **returns its file
name** (a string), which is what makes `Close` convenient as the last line of a
round trip. `Close["file"]` / `Close[File["file"]]` closes the stream open for that
name — the most recently opened, if several share it.

Closing frees the registry slot and, for an output stream, flushes and closes the
underlying file. Closing a stream that is not open returns `$Failed` with a
`Close::stream` message. Streams left open are closed automatically at program
exit, so a forgotten `Close` does not leak — but an output file is only guaranteed
flushed to disk once it is closed.
