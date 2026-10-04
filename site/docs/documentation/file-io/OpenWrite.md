# OpenWrite

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`OpenWrite["file"]`**

opens a file for writing (truncating it) and returns an OutputStream object; also accepts File\["file"\]. Returns $Failed on failure.

## Examples (15)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (9)

OpenWrite truncates the file

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_io_open.txt"]
Out[1]= OutputStream["/tmp/mathilda_io_open.txt", 1]
```

```mathematica
In[2]:= Write[str, x^2]

In[3]:= Close[str]
Out[3]= "/tmp/mathilda_io_open.txt"
```

OpenAppend keeps the prior contents

```mathematica
In[4]:= app = OpenAppend["/tmp/mathilda_io_open.txt"]
Out[4]= OutputStream["/tmp/mathilda_io_open.txt", 2]
```

```mathematica
In[5]:= Write[app, y^2]

In[6]:= Close[app]
Out[6]= "/tmp/mathilda_io_open.txt"

In[7]:= ins = OpenRead["/tmp/mathilda_io_open.txt"]
Out[7]= InputStream["/tmp/mathilda_io_open.txt", 3]
```

Both lines survive the append

```mathematica
In[8]:= ReadList[ins]
Out[8]= {x^2, y^2}
```

```mathematica
In[9]:= Close[ins]
Out[9]= "/tmp/mathilda_io_open.txt"
```

### Applications (6)

Open for writing (object bound to str)

```mathematica
In[10]:= str = OpenWrite["/tmp/mathilda_ow.txt"];
```

Write two lines, then close

```mathematica
In[11]:= Write[str, x + y]; Write[str, x^2]; Close[str];
```

Read them back

```mathematica
In[12]:= ReadList["/tmp/mathilda_ow.txt"]
Out[12]= {x + y, x^2}
```

```mathematica
In[13]:= str = OpenWrite["/tmp/mathilda_ow2.txt"]; Write[str, 111]; Close[str];
```

OpenWrite truncates: 111 is gone

```mathematica
In[14]:= str = OpenWrite["/tmp/mathilda_ow2.txt"]; Write[str, 222]; Close[str];
```

```mathematica
In[15]:= ReadList["/tmp/mathilda_ow2.txt"]
Out[15]= {222}
```

## Implementation notes

**Algorithm.** `builtin_openwrite` calls `open_common(res, output=1, append=0, "OpenWrite")`,
which resolves the name (`stream_filename_arg` accepts `"file"` or `File["file"]`) and calls
`stream_open_output(name, append=0)` — `fopen(name, "wb")`, which **truncates** any existing
file. A fresh registry slot records `{name, fp, is_output = 1}` and the result is the inert
object `OutputStream["file", id]`. A file that cannot be opened prints `OpenWrite::noopen` and
returns `$Failed`.

**Data structures — the stream layer.** Shares the process-global `Stream` registry described
under `OpenRead`. An output stream holds an open `FILE*` rather than a buffer;
`id` in `OutputStream["file", id]` indexes the slot so `Write` /
`StreamPosition` address it across calls. The `atexit` hook `fclose`s any
stream still open at exit, so an output file is flushed even if `Close` is
forgotten — though `Write`/`WriteString` `fflush` after every
call, so a round trip with `Read` sees the bytes immediately.

**Complexity / limits.** `O(1)`; the open is a single `fopen`. The truncation is the
distinction from `OpenAppend` (`"ab"`). Pure ANSI C99 (`fopen`), no POSIX
guards. `ATTR_PROTECTED`.

- `Protected`. Return `$Failed` (with an `Open*::noopen` diagnostic) if the file cannot be opened.
- The integer in the returned object is an internal handle into the stream registry; the object is inert and prints as itself.
- The registry is freed on `Close` or, for anything still open, at program exit (no leaks).

**Attributes:** `Protected`.

## References

**See also:** [OpenRead](../../file-io/OpenRead/), [OpenAppend](../../file-io/OpenAppend/), [Close](../../file-io/Close/)

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`OpenWrite["file"]` opens a file for writing and returns an
`OutputStream["file", n]` object; the integer `n` is an internal handle into the
stream registry. It **truncates** any existing file, so the second block above ends
with just `222`. To keep prior contents and add to the end, use
`OpenAppend` instead.

A named file that is not already open is also auto-opened (truncating) by
`Write`/`WriteString`, so an explicit `OpenWrite` is
only needed to hold the stream across several writes or to interleave with
`StreamPosition`. Finish with `Close`, which returns the file name; any
stream still open is closed at program exit, so nothing leaks.
