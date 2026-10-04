# Streams

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Streams[]`**

gives a list of all currently open streams.

**`Streams["file"] lists only the open streams for the named file.`**

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (6)

```mathematica
In[1]:= a = OpenWrite["/tmp/mathilda_io_a.txt"]
Out[1]= OutputStream["/tmp/mathilda_io_a.txt", 1]

In[2]:= b = OpenWrite["/tmp/mathilda_io_b.txt"]
Out[2]= OutputStream["/tmp/mathilda_io_b.txt", 2]
```

All open streams

```mathematica
In[3]:= Streams[]
Out[3]= {OutputStream["/tmp/mathilda_io_a.txt", 1], OutputStream["/tmp/mathilda_io_b.txt", 2]}
```

Just the streams for one file

```mathematica
In[4]:= Streams["/tmp/mathilda_io_a.txt"]
Out[4]= {OutputStream["/tmp/mathilda_io_a.txt", 1]}
```

```mathematica
In[5]:= Close[a]
Out[5]= "/tmp/mathilda_io_a.txt"

In[6]:= Close[b]
Out[6]= "/tmp/mathilda_io_b.txt"
```

### Applications (2)

The open stream is listed

```mathematica
In[7]:= str = OpenWrite["/tmp/mathilda_streams.txt"]; MemberQ[Streams[], str]
Out[7]= True
```

After Close, none for that file

```mathematica
In[8]:= Close[str]; Streams["/tmp/mathilda_streams.txt"]
Out[8]= {}
```

## Implementation notes

**Algorithm.** `builtin_streams` walks the process-global stream registry and builds a `List`
of the live stream objects. With one argument it reads a file-name filter
(`stream_filename_arg`) and keeps only slots whose `name` matches it; with none it keeps every
open slot. For each surviving slot it emits `stream_make_object` (`InputStream["file", id]` or
`OutputStream["file", id]` by the slot's `is_output` flag) into a `malloc`'d array that is
grown as needed, then wraps it in `List[...]`.

**Data structures — the stream layer.** Reads the same process-global `Stream` registry
described under `OpenRead`; the result is a snapshot of the slots whose
`id != 0`, so it shrinks as `Close` frees slots and is empty at program start and
once everything is closed. Each entry is the inert object itself, so it can be passed straight
back to `Read`/`Write`/`Close`.

**Complexity / limits.** `O(n)` over the registry, `n` the number of open slots. The id inside
each object is an internal handle whose value depends on how many streams have been opened, so
code should test membership rather than match a printed form. `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## References

**See also:** [InputStream](../../other-advanced/InputStream/), [OutputStream](../../other-advanced/OutputStream/)

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`Streams[]` returns the list of currently open `InputStream`/`OutputStream`
objects; `Streams["file"]` lists only those open for the named file. Each entry is
the same inert object `OpenRead`/`OpenWrite` handed
back, so it can be passed straight to `Read`, `Write` or
`Close`.

The list reflects the live registry, so it shrinks as streams are closed and is
empty once everything is closed (and at program start). The integer in each object
is an internal handle whose exact value depends on how many streams have been
opened in the session, so membership (`MemberQ[Streams[], str]`) is the robust way
to test for a stream rather than matching its printed form.
