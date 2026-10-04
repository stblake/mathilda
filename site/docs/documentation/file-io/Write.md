# Write

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Write[stream, expr1, expr2, ...]`**

writes the expressions to an output stream in input form, followed by a newline. The stream may be an OutputStream, a "file", or File\["file"\]; a named file that is not already open is opened for writing and left open.

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (8)

```mathematica
In[1]:= str = OpenWrite["/tmp/mathilda_io_w.txt"]
Out[1]= OutputStream["/tmp/mathilda_io_w.txt", 1]
```

Write evaluates, prints input form, adds a newline

```mathematica
In[2]:= Write[str, 1 + 1]
```

```mathematica
In[3]:= Write[str, a + b]
```

WriteString writes raw text, no added newline

```mathematica
In[4]:= WriteString[str, "done\n"]
```

```mathematica
In[5]:= Close[str]
Out[5]= "/tmp/mathilda_io_w.txt"
```

1 + 1 was evaluated to 2 before writing

```mathematica
In[6]:= ReadList["/tmp/mathilda_io_w.txt"]
Out[6]= {2, a + b, done}
```

Hold writes the unevaluated form

```mathematica
In[7]:= str = OpenWrite["/tmp/mathilda_io_wh.txt"]; Write[str, Hold[2 + 2]]; Close[str];
```

```mathematica
In[8]:= ReadList["/tmp/mathilda_io_wh.txt", String]
Out[8]= {"Hold[2 + 2]"}
```

### Applications (5)

```mathematica
In[9]:= str = OpenWrite["/tmp/mathilda_write.txt"];
```

Each Write evaluates, then adds a newline

```mathematica
In[10]:= Write[str, 1 + 1]; Write[str, a + b]; Close[str];
```

1 + 1 was evaluated to 2 before writing

```mathematica
In[11]:= ReadList["/tmp/mathilda_write.txt"]
Out[11]= {2, a + b}
```

Hold writes the unevaluated form

```mathematica
In[12]:= str = OpenWrite["/tmp/mathilda_write2.txt"]; Write[str, Hold[2 + 2]]; Close[str];
```

```mathematica
In[13]:= ReadList["/tmp/mathilda_write2.txt", String]
Out[13]= {"Hold[2 + 2]"}
```

## Implementation notes

**Algorithm.** `builtin_write` resolves its first argument with `resolve_output`: an
`OutputStream` object, or a `"file"`/`File["file"]` name that is **auto-opened** (truncating)
and left open; an `InputStream` or an unopenable name fails with `$Failed` (`Write::openx` /
`Write::noopen`). For each remaining argument it renders the expression with `expr_to_string`
(re-readable input form) and `fputs` it, then writes a single `'\n'` and `fflush`es. Returns
`Null`. Because its arguments are ordinary (evaluated) arguments, `1 + 1` is written as `2`;
`Hold[...]` is the way to write an unevaluated form.

**Data structures — the stream layer.** Writes to the `FILE*` of an output slot in the
process-global `Stream` registry described under `OpenRead`. The `fflush` after
each call is what makes the bytes immediately visible to `Read`/`ReadList`,
so a write-then-read round trip works without an intervening `Close`.

**Complexity / limits.** `O(total output length)`. The trailing newline per call is what makes
the output parse back one object per `Read`. Differs from `WriteString`,
which writes strings raw with no quoting and no added newline. Pure ANSI C99
(`fputs`/`fputc`/`fflush`). `ATTR_PROTECTED`.

- `Protected`. Return `$Failed` if the file cannot be opened.
- `Write` evaluates its expression arguments before writing (use `Hold[...]` to write an unevaluated form), and output is flushed after each call so it round-trips with `Read`/`ReadList`.

**Attributes:** `Protected`.

## References

**See also:** [WriteString](../../file-io/WriteString/), [OutputStream](../../other-advanced/OutputStream/), [Read](../../file-io/Read/), [ReadList](../../file-io/ReadList/)

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`Write[stream, e1, e2, ...]` writes each expression in re-readable input form,
followed by a single newline, and flushes — so the output round-trips through
`Read`/`ReadList`. Expression arguments are **evaluated**
first (`1 + 1` becomes `2`); wrap in `Hold[...]` to write the unevaluated form.

`stream` may be an `OutputStream`, a `"file"` string, or `File["file"]`; a named
file that is not already open is auto-opened (truncating) and left open. Writing
raw text with no quoting or added newline is `WriteString`'s job.
`Write` returns `Null`.
