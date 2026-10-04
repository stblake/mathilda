# OpenRead

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`OpenRead["file"]`**

opens a file for reading and returns an InputStream object; also accepts File\["file"\]. Returns $Failed if the file cannot be opened.

## Examples (14)

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

### Applications (5)

Make a file

```mathematica
In[10]:= str = OpenWrite["/tmp/mathilda_or.txt"]; Write[str, 10]; Write[str, 20]; Close[str];
```

Open it for reading

```mathematica
In[11]:= ins = OpenRead["/tmp/mathilda_or.txt"];
```

First object

```mathematica
In[12]:= Read[ins, Number]
Out[12]= 10
```

The current point has advanced

```mathematica
In[13]:= Read[ins, Number]
Out[13]= 20
```

```mathematica
In[14]:= Close[ins]
Out[14]= "/tmp/mathilda_or.txt"
```

## Implementation notes

**Algorithm.** `builtin_openread` calls `open_common(res, output=0, append=0, "OpenRead")`,
which pulls the name out of a `"file"` string or `File["file"]` (`stream_filename_arg`) and
calls `stream_open_input`. That **slurps the whole file** into a `malloc`'d buffer
(`slurp_file`: `fopen("rb")`, `fseek`/`ftell` for the size, one `fread`, NUL-terminate),
claims a registry slot, and records `{name, buf, len, pos = 0}`. The handed-back value is the
inert object `InputStream["file", id]` built by `stream_make_object`, where `id` is the slot's
monotonic handle. A file that cannot be opened prints `OpenRead::noopen` and returns `$Failed`.

**Data structures — the stream layer.** `src/io/streams.c` keeps one **process-global stream
registry**: a dynamically grown array of `Stream` slots (`alloc_slot` doubles it; a slot with
`id == 0` is free), with `g_next_id` a monotonic handle allocator so a reused slot never
reuses a live id. An input stream is a resident buffer plus a moving **current point** `pos`;
the integer in `InputStream["file", id]` indexes the registry across calls, which is how
successive `Read` / `StreamPosition` calls address the same
stream. A `streams_shutdown` `atexit` hook frees every slot's `name`/`buf` and `fclose`s any
open `fp`, so nothing leaks under valgrind.

**Complexity / limits.** `O(file size)` time and memory for the initial slurp (the entire file
is resident, so a very large file is fully loaded), then `O(1)` per seek. Pure ANSI C99
(`fopen`/`fread`/`fseek`/`ftell`), so no POSIX feature-test guards are needed. `ATTR_PROTECTED`.

- `Protected`. Return `$Failed` (with an `Open*::noopen` diagnostic) if the file cannot be opened.
- The integer in the returned object is an internal handle into the stream registry; the object is inert and prints as itself.
- The registry is freed on `Close` or, for anything still open, at program exit (no leaks).

**Attributes:** `Protected`.

## References

**See also:** [OpenWrite](../../file-io/OpenWrite/), [OpenAppend](../../file-io/OpenAppend/), [Close](../../file-io/Close/)

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`OpenRead["file"]` (or `OpenRead[File["file"]]`) opens a file for reading and
returns an `InputStream["file", n]` object. The whole file is slurped into a buffer
with a moving **current point**, so successive `Read` calls on the same
stream return successive objects, and `StreamPosition` /
`SetStreamPosition` query and move that point.

`$Failed` (with an `OpenRead::noopen` diagnostic) is returned when the file cannot
be opened. A named file handed directly to `Read`/`ReadList` is auto-opened, so an
explicit `OpenRead` is for when you want to hold the stream and advance it yourself.
Finish with `Close`, which returns the file name.
