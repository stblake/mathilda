# WriteString

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`WriteString[stream, str1, str2, ...]`**

writes the strings to an output stream with no added quotes or newline. Non-string arguments are written in input form.

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
In[9]:= str = OpenWrite["/tmp/mathilda_ws.txt"];
```

Raw: no quotes, no separators, no newline

```mathematica
In[10]:= WriteString[str, "a", "b", "c"]; Close[str];
```

The whole line is "abc"

```mathematica
In[11]:= ReadList["/tmp/mathilda_ws.txt", String]
Out[11]= {"abc"}
```

Embed newlines yourself

```mathematica
In[12]:= str = OpenWrite["/tmp/mathilda_ws2.txt"]; WriteString[str, "1 2 3\n4 5 6\n"]; Close[str];
```

```mathematica
In[13]:= ReadList["/tmp/mathilda_ws2.txt", Number]
Out[13]= {1, 2, 3, 4, 5, 6}
```

## Implementation notes

**Algorithm.** `builtin_writestring` resolves its first argument with `resolve_output` exactly
as `Write` does (an `OutputStream`, or a name auto-opened truncating and left open;
`$Failed` otherwise). It then writes each remaining argument **verbatim**: a string is `fputs`'d
raw, with no surrounding quotes, and a non-string argument is rendered with `expr_to_string`
first. It adds **no separators between arguments and no trailing newline**, then `fflush`es and
returns `Null`.

**Data structures — the stream layer.** Writes to the `FILE*` of an output slot in the
process-global `Stream` registry described under `OpenRead`. It is the exact-byte
counterpart to `Write`: `Write` prints input form and appends `'\n'`, whereas
`WriteString` emits precisely the bytes given — so a caller lays out CSV lines, headers, or any
newline structure explicitly (e.g. `WriteString[s, "1 2 3\n4 5 6\n"]`).

**Complexity / limits.** `O(total output length)`; `fflush` per call. Pure ANSI C99
(`fputs`/`fflush`). `ATTR_PROTECTED`.

- `Protected`. Return `$Failed` if the file cannot be opened.
- `Write` evaluates its expression arguments before writing (use `Hold[...]` to write an unevaluated form), and output is flushed after each call so it round-trips with `Read`/`ReadList`.

**Attributes:** `Protected`.

## References

**See also:** [Write](../../file-io/Write/), [OutputStream](../../other-advanced/OutputStream/), [Read](../../file-io/Read/), [ReadList](../../file-io/ReadList/)

- Source: [`src/io/streams.c`](https://github.com/stblake/mathilda/blob/main/src/io/streams.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_streams.c`](https://github.com/stblake/mathilda/blob/main/tests/test_streams.c)

## Notes & additional examples

### Notes

`WriteString[stream, s1, s2, ...]` writes the strings **verbatim** — no surrounding
quotes, no separators between arguments, and no trailing newline. It is the tool
for laying out an exact byte stream (a CSV line, a header, a line you terminate
with an explicit `"\n"`), in contrast to `Write`, which prints input
form and appends a newline.

A non-string argument is written in input form. Output is flushed after each call.
`stream` may be an `OutputStream`, a `"file"`, or `File["file"]`; a named file that
is not already open is auto-opened (truncating) and left open. `WriteString`
returns `Null`.
