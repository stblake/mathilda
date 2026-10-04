# FileSize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FileSize["name"]`**

gives the number of bytes in the file with the specified name.

<details>
<summary>Notes</summary>

In FileSize\["name"\], name is interpreted relative to your current directory. FileSize does not search $Path. FileSize follows symbolic links, reporting the size of the target file. FileSize gives the size as an integer count of bytes, not a Quantity. FileSize returns $Failed and prints a message if the file cannot be found.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

Writes "2\n3\n5\n7\n11\n" = 11 bytes

```mathematica
In[1]:= Put[2, 3, 5, 7, 11, "/tmp/mathilda_io_fs.txt"]
```

```mathematica
In[2]:= FileSize["/tmp/mathilda_io_fs.txt"]
Out[2]= 11
```

A plain Integer, not a Quantity

```mathematica
In[3]:= Head[FileSize["/tmp/mathilda_io_fs.txt"]]
Out[3]= Integer
```

Missing file -> $Failed

```mathematica
In[4]:= Quiet[FileSize["/tmp/mathilda_no_such_file_xyz"]]
Out[4]= $Failed
```

### Applications (4)

Writes "2\n3\n5\n7\n11\n" = 11 bytes

```mathematica
In[5]:= Put[2, 3, 5, 7, 11, "/tmp/mathilda_fs.txt"];
```

```mathematica
In[6]:= FileSize["/tmp/mathilda_fs.txt"]
Out[6]= 11
```

A plain Integer, not a Quantity

```mathematica
In[7]:= Put[2, 3, 5, 7, 11, "/tmp/mathilda_fs.txt"]; Head[FileSize["/tmp/mathilda_fs.txt"]]
Out[7]= Integer
```

Missing file -> $Failed

```mathematica
In[8]:= Quiet[FileSize["/tmp/mathilda_no_such_file_xyz"]]
Out[8]= $Failed
```

## Implementation notes

**Algorithm.** `builtin_filesize` accepts exactly one `EXPR_STRING` argument (anything else —
wrong arity, a symbol, a non-string atom — returns `NULL`, leaving the call unevaluated so
symbolic arguments flow through). It calls `stat()` on the path and, on success, returns
`st.st_size` as an `Integer` (`expr_new_integer((int64_t)st.st_size)`) — a plain byte count,
not a `Quantity`. `stat` (not `lstat`) follows symbolic links, so the size reported is the
target file's.

When the path cannot be `stat`'d (most often because nothing is there) it emits
`FileSize::nffil` through `fs_msg` and returns `$Failed`. `fs_msg` is the local diagnostic
funnel: it calls `mth_msg_note_fired()` (so `Check[]` sees the message) and returns early under
`mth_msg_suppressed()` (so `Quiet[]` silences it), then `vfprintf`s to `stderr` — the same
`Quiet`/`Check`-aware pattern as `dt_msg` in `src/datetime.c`.

**Data structures.** None beyond a `struct stat`; the path is read directly out of the argument
string. `_POSIX_C_SOURCE 200809L` is defined before any include so `stat` is visible under
glibc's `-std=c99` (SPEC.md §10).

**Complexity / limits.** `O(1)` — one `stat` syscall. The path is interpreted relative to the
current working directory; `$Path` is not searched. `ATTR_PROTECTED`.

- `Protected`.
- `"name"` is interpreted relative to the current working directory. `$Path` is not searched.
- Implemented with `stat()`, so symbolic links are followed and the target file's size is reported.
- Returns `$Failed` and prints a `FileSize::nffil` message when the file cannot be found. The message respects `Quiet[]` and is visible to `Check[]`.
- Leaves the call unevaluated when given the wrong arity, a symbolic argument, or any non-string atom.

**Attributes:** `Protected`.

## References

- Source: [`src/files.c`](https://github.com/stblake/mathilda/blob/main/src/files.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_files.c`](https://github.com/stblake/mathilda/blob/main/tests/test_files.c)

## Notes & additional examples

### Notes

`FileSize["name"]` gives the file's size in **bytes** as a plain `Integer` (not a
`Quantity`). The path is relative to the current working directory; `$Path` is not
searched. It is implemented with `stat`, so a symbolic link is followed and the
target file's size is reported.

A path that cannot be `stat`'d (usually because nothing is there) returns `$Failed`
with a `FileSize::nffil` message, which respects `Quiet[]` and is visible to
`Check[]`. A symbolic or non-string argument leaves the call unevaluated.
