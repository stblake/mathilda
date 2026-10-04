# File

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`File["name"]`**

is a symbolic wrapper for a file name, accepted wherever Read, ReadList, OpenRead, OpenWrite, OpenAppend, and Close take a file.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

An inert wrapper head, kept unevaluated

```mathematica
In[1]:= Head[File["data.txt"]]
Out[1]= File
```

A one-argument File[...] holding the path string

```mathematica
In[2]:= FullForm[File["data.txt"]]
Out[2]= File["data.txt"]
```

```mathematica
In[3]:= str = OpenWrite["/tmp/mathilda_word.txt"]; WriteString[str, "alpha beta gamma\n"]; Close[str];
```

A File wrapper works like a bare path

```mathematica
In[4]:= ReadList[File["/tmp/mathilda_word.txt"], Word]
Out[4]= {"alpha", "beta", "gamma"}
```

## Implementation notes

**Definition.** `File["name"]` is a **symbolic wrapper** for a file name, accepted wherever
`Read`, `ReadList`, `OpenRead`, `OpenWrite`, `OpenAppend` and `Close` take a file. It is an
inert head, not a function: `streams_init` in `src/io/streams.c` marks `File` `Protected`
(the docstring is set centrally in `src/info.c`), and there is no builtin that rewrites it —
`File["x"]` evaluates to itself.

**Representation.** `File["name"]` is an ordinary `EXPR_FUNCTION` with head `SYM_File` and
a single string argument, so `Head[File["data.txt"]]` is `File` and
`FullForm[File["data.txt"]]` is `File["data.txt"]`. The I/O layer never stores a `File`
object itself; the helper `stream_filename_arg` (`src/io/streams.c`) unwraps it, accepting
either a bare string or a one-argument `File[...]` whose argument is a string, and returns
the underlying path.

**Usage & limits.** `File` lets a path be passed as an explicit, typed object rather than a
raw string — `ReadList[File["data.txt"], Word]` behaves exactly like
`ReadList["data.txt", Word]`. The argument must be a single string; other forms are not
unwrapped. It is a name wrapper only: it does not open, test or resolve the file, and holds
no handle (an open stream is an `InputStream`/`OutputStream` object instead).

**Attributes:** `Protected`.

## References

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`File["name"]` is a symbolic wrapper for a file name, accepted wherever `Read`, `ReadList`,
`OpenRead`, `OpenWrite`, `OpenAppend` and `Close` take a file. It is an inert, `Protected`
head: `File["x"]` evaluates to itself, and `Head[File["x"]]` is `File`.

The I/O layer unwraps it on the way in — `stream_filename_arg` accepts either a bare string
or a one-argument `File[...]` and takes the path from inside — so
`ReadList[File["data.txt"], Word]` is exactly `ReadList["data.txt", Word]`. It is a name
wrapper only: it holds no handle and does not open or test the file (an open stream is an
`InputStream`/`OutputStream` object).
