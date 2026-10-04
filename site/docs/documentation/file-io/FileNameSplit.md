# FileNameSplit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FileNameSplit["name"]`**

splits a file name into a list of parts.

**`FileNameSplit[..., OperatingSystem->"os"] uses the conventions of the specified operating system; possible choices are "Windows", "MacOSX", and "Unix".`**

<details>
<summary>Notes</summary>

FileNameSplit by default uses pathname separators and other conventions suitable for your operating system. Absolute file names that begin with a pathname separator yield a list of parts that starts with "". Under Windows, the drive or share name is treated as the first part of the file name. FileNameSplit just operates on names of files; it does not actually search for the file specified.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

Leading "" marks an absolute path

```mathematica
In[1]:= FileNameSplit["/home/user/data.csv"]
Out[1]= {"", "home", "user", "data.csv"}
```

Duplicate and trailing separators are dropped

```mathematica
In[2]:= FileNameSplit["a//b/c/"]
Out[2]= {"a", "b", "c"}
```

Join inverts Split

```mathematica
In[3]:= FileNameJoin[FileNameSplit["/usr/local/bin"]]
Out[3]= "/usr/local/bin"
```

A drive is an ordinary first part

```mathematica
In[4]:= FileNameSplit["C:\\path\\file.txt", OperatingSystem -> "Windows"]
Out[4]= {"C:", "path", "file.txt"}
```

### Applications (4)

Leading "" marks an absolute path

```mathematica
In[5]:= FileNameSplit["/home/user/data.csv"]
Out[5]= {"", "home", "user", "data.csv"}
```

Duplicate and trailing separators are dropped

```mathematica
In[6]:= FileNameSplit["a//b/c/"]
Out[6]= {"a", "b", "c"}
```

Join inverts Split

```mathematica
In[7]:= FileNameJoin[FileNameSplit["/usr/local/bin"]]
Out[7]= "/usr/local/bin"
```

A drive is an ordinary first part

```mathematica
In[8]:= FileNameSplit["C:\\path\\file.txt", OperatingSystem -> "Windows"]
Out[8]= {"C:", "path", "file.txt"}
```

## Implementation notes

**Algorithm.** `builtin_filenamesplit` is the structural inverse of
`FileNameJoin` and, like it, is a **pure string operation** that never
touches the filesystem. It decodes the same trailing `OperatingSystem -> "..."` option into a
`windows` flag (sharing `fnj_is_sep` and the absolute/UNC rules), requires a single string
`spec`, and calls `fns_split_build`.

`fns_split_build` handles the leading context first: under Windows a `\\host\share` UNC prefix
is captured as one part; otherwise a leading separator (an **absolute path**) emits a leading
`""` part. It then consumes the rest as maximal non-separator runs, dropping empty runs from
trailing and duplicated separators — so `"a//b/c/"` splits to `{"a","b","c"}`. The parts are
assembled into a `List[...]`; a non-string argument or an unknown OS leaves the call
unevaluated, and `FileNameSplit[]` prints `FileNameSplit::argx`.

**Data structures.** A `char**` of freshly-`malloc`'d part strings (`fns_dup` copies each
substring; an upper bound of `strlen(path)+2` parts is allocated once), converted to
`expr_new_string` elements that `expr_new_function` adopts into the `List`. A Windows drive like
`C:` contains no separator and so falls out naturally as an ordinary first part.

**Complexity / limits.** `O(path length)`, single pass. `FileNameJoin[FileNameSplit[name]]`
reconstructs a canonicalized `name`. `ATTR_PROTECTED`.

- `Protected`.
- Pure string operation — does not touch the filesystem.
- A leading pathname separator marks an absolute path and yields a leading `""` part; trailing and duplicate separators are dropped (`"a//b/"` → `{"a", "b"}`).
- `"Windows"` treats a leading `\\server\share` UNC prefix as a single part and a drive like `C:` as an ordinary first part; `"MacOSX"`/`"Unix"` split on `/`. The default is the host operating system's separator.
- `Options[FileNameSplit]` reports the `OperatingSystem` default.
- `FileNameJoin[FileNameSplit[name]]` reconstructs a canonicalized `name`.
- `FileNameSplit[]` prints `FileNameSplit::argx` and stays unevaluated; a non-string argument or an unknown OS leaves the call unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [FileNameJoin](../../file-io/FileNameJoin/)

- Source: [`src/files.c`](https://github.com/stblake/mathilda/blob/main/src/files.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_files.c`](https://github.com/stblake/mathilda/blob/main/tests/test_files.c)

## Notes & additional examples

### Notes

`FileNameSplit["name"]` is the structural inverse of
`FileNameJoin`: it returns the list of path components. It is a
**pure string operation** and never touches the filesystem. A leading separator
makes the path absolute and yields a leading `""` part; trailing and duplicate
separators are dropped.

The separator defaults to the host operating system's; `OperatingSystem ->
"Windows" | "MacOSX" | "Unix"` selects it. On `"Windows"` a leading UNC
`\\server\share` prefix is kept as a single part and a drive like `C:` falls out
as an ordinary first part. A non-string argument or an unknown OS leaves the call
unevaluated.
