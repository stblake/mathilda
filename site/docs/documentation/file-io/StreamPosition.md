# StreamPosition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`StreamPosition[stream]`**

gives the position of the current point in an open stream, as an integer byte offset.

## Examples (13)

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

### Applications (5)

Writes "100\n200\n300\n"

```mathematica
In[9]:= Put[100, 200, 300, "/tmp/mathilda_sp.txt"];
```

```mathematica
In[10]:= ins = OpenRead["/tmp/mathilda_sp.txt"];
```

100

```mathematica
In[11]:= Read[ins, Number]
Out[11]= 100
```

Byte offset after "100\n" is 4

```mathematica
In[12]:= StreamPosition[ins]
Out[12]= 4
```

```mathematica
In[13]:= Close[ins]
Out[13]= "/tmp/mathilda_sp.txt"
```

## Implementation notes

**Algorithm.** `builtin_streamposition` resolves its argument with `resolve_open` — a stream
object (`stream_handle`) or an already-open file name (`stream_find_by_name`), with **no
auto-open**: a stream that is not already open prints `StreamPosition::stream` and the call
returns `$Failed`. The position is read as a byte offset: for an output stream it is
`ftell(s->fp)`, for an input stream it is the buffer cursor `(long)s->pos`. A negative `ftell`
is clamped to `0`, and the offset is returned as an `Integer`.

**Data structures — the stream layer.** Reads the current point of a slot in the process-global
`Stream` registry described under `OpenRead`: an input stream's `pos` into its
resident buffer, or an output stream's file offset. It is the query half of the seek pair whose
setter is `SetStreamPosition`; reading through `Read`
advances this same `pos`.

**Complexity / limits.** `O(n)` over the registry to locate the slot, then `O(1)`. The offset
is in **bytes**, not characters or objects. Pure ANSI C99 (`ftell`). `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## References

**See also:** [SetStreamPosition](../../file-io/SetStreamPosition/)

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`StreamPosition[stream]` returns the stream's current point as an integer **byte
offset** from the start of the file. For an input stream it is the position in the
slurped buffer; for an output stream it is `ftell` of the underlying file. It pairs
with `SetStreamPosition`, which moves the point, so the two
together give random access within a file.

Reading advances the offset, which is why the value above is `4` after one
`Read` of `"100\n"`. `$Failed` is returned for a stream that is not open.
