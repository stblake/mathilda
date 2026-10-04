# SetStreamPosition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SetStreamPosition[stream, n]`**

sets the current point of a stream to byte offset n and returns the new position. SetStreamPosition\[stream, Infinity\] moves to the end.

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (8)

Writes "100\n200\n300\n"

```mathematica
In[1]:= Put[100, 200, 300, "/tmp/mathilda_io_sp.txt"]
```

```mathematica
In[2]:= str = OpenRead["/tmp/mathilda_io_sp.txt"]
Out[2]= InputStream["/tmp/mathilda_io_sp.txt", 1]

In[3]:= Read[str, Number]
Out[3]= 100
```

Byte offset after "100\n"

```mathematica
In[4]:= StreamPosition[str]
Out[4]= 4
```

Rewind to the start

```mathematica
In[5]:= SetStreamPosition[str, 0]
Out[5]= 0
```

The same first number again

```mathematica
In[6]:= Read[str, Number]
Out[6]= 100
```

Jump to end of file: the byte length

```mathematica
In[7]:= SetStreamPosition[str, Infinity]
Out[7]= 12
```

```mathematica
In[8]:= Close[str]
Out[8]= "/tmp/mathilda_io_sp.txt"
```

### Applications (7)

```mathematica
In[9]:= Put[100, 200, 300, "/tmp/mathilda_ssp.txt"];

In[10]:= ins = OpenRead["/tmp/mathilda_ssp.txt"];
```

100

```mathematica
In[11]:= Read[ins, Number]
Out[11]= 100
```

Rewind to the start; returns the new position

```mathematica
In[12]:= SetStreamPosition[ins, 0]
Out[12]= 0
```

100 again

```mathematica
In[13]:= Read[ins, Number]
Out[13]= 100
```

Jump to end of file: the byte length

```mathematica
In[14]:= SetStreamPosition[ins, Infinity]
Out[14]= 12
```

```mathematica
In[15]:= Close[ins]
Out[15]= "/tmp/mathilda_ssp.txt"
```

## Implementation notes

**Algorithm.** `builtin_setstreamposition` resolves an already-open stream with `resolve_open`
(no auto-open; `$Failed` + `SetStreamPosition::stream` otherwise), then interprets its second
argument. `Infinity` means the end of the stream; otherwise the target must be an `Integer`,
and a negative value is clamped to `0`. For an **output** stream it `fseek`s (`SEEK_END` for
`Infinity`, else `SEEK_SET` to the offset) and returns `ftell`. For an **input** stream it sets
the buffer cursor `s->pos`, clamping a past-the-end target to the buffer length `s->len`
(`Infinity` → `len`). Either way it returns the resulting byte offset as an `Integer`.

**Data structures — the stream layer.** Moves the current point of a slot in the process-global
`Stream` registry described under `OpenRead`. It is the setter half of the seek
pair whose getter is `StreamPosition`; the clamping keeps the point within
`[0, len]` for an input stream, so a rewind-and-reread or skip-ahead is always in bounds.

**Complexity / limits.** `O(n)` over the registry to locate the slot, then `O(1)`. Offsets are
in **bytes**. A non-integer, non-`Infinity` target leaves the call unevaluated. Pure ANSI C99
(`fseek`/`ftell`). `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## References

**See also:** [StreamPosition](../../file-io/StreamPosition/)

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`SetStreamPosition[stream, n]` moves the current point to byte offset `n` and
returns the new position; `SetStreamPosition[stream, Infinity]` moves to the end of
the stream (so it reports the byte length). A negative `n` clamps to `0` and an `n`
past the end clamps to the length, so the position stays within the file.

Combined with `StreamPosition` this gives seek-style random
access: rewind and re-read, or skip ahead. `$Failed` is returned for a stream that
is not open.
