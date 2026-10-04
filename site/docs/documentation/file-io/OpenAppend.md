# OpenAppend

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`OpenAppend["file"]`**

opens a file for appending and returns an OutputStream object; also accepts File\["file"\]. Returns $Failed on failure.

## Examples (13)

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

### Applications (4)

File now holds 1

```mathematica
In[10]:= str = OpenWrite["/tmp/mathilda_oa.txt"]; Write[str, 1]; Close[str];
```

Open at the end, keeping the 1

```mathematica
In[11]:= app = OpenAppend["/tmp/mathilda_oa.txt"];
```

```mathematica
In[12]:= Write[app, 2]; Write[app, 3]; Close[app];
```

All three lines survive

```mathematica
In[13]:= ReadList["/tmp/mathilda_oa.txt"]
Out[13]= {1, 2, 3}
```

## Implementation notes

**Algorithm.** `builtin_openappend` calls `open_common(res, output=1, append=1, "OpenAppend")`,
which resolves the name and calls `stream_open_output(name, append=1)` — `fopen(name, "ab")`.
Unlike `OpenWrite`'s `"wb"`, the `"ab"` mode **preserves** the existing
contents and positions writes at the end, and creates the file if it does not exist. A slot
records `{name, fp, is_output = 1}` and the result is `OutputStream["file", id]`. A file that
cannot be opened prints `OpenAppend::noopen` and returns `$Failed`.

**Data structures — the stream layer.** Shares the process-global `Stream` registry described
under `OpenRead`; the only difference from `OpenWrite` is the `fopen` mode, so
the returned object, the id handle, and the `Close`/`atexit` lifecycle are identical.

**Complexity / limits.** `O(1)`. Pure ANSI C99 (`fopen`), no POSIX guards. `ATTR_PROTECTED`.

- `Protected`. Return `$Failed` (with an `Open*::noopen` diagnostic) if the file cannot be opened.
- The integer in the returned object is an internal handle into the stream registry; the object is inert and prints as itself.
- The registry is freed on `Close` or, for anything still open, at program exit (no leaks).

**Attributes:** `Protected`.

## References

**See also:** [OpenRead](../../file-io/OpenRead/), [OpenWrite](../../file-io/OpenWrite/), [Close](../../file-io/Close/)

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`OpenAppend["file"]` opens a file for writing **at its end** and returns an
`OutputStream["file", n]` object. Unlike `OpenWrite`, which
truncates, it preserves the existing contents and adds after them — the round trip
above ends with `{1, 2, 3}`. If the file does not exist it is created empty, so
`OpenAppend` is safe on a first run.

`$Failed` (with an `OpenAppend::noopen` diagnostic) is returned when the file
cannot be opened. Finish with `Close`, which returns the file name.
