# FileNameJoin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FileNameJoin[{"name1", "name2", ...}]`**

joins the namei into a file name suitable for your current operating system.

**`FileNameJoin[{"", "name1", ...}] gives an absolute file path beginning with a pathname separator.`**

**`FileNameJoin["name"] canonicalizes name, making pathname separators appropriate for your operating system.`**

**`FileNameJoin[..., OperatingSystem->"os"] yields a file name in the format for the specified operating system; possible choices are "Windows", "MacOSX", and "Unix".`**

<details>
<summary>Notes</summary>

The namei can be individual names or file paths containing pathname separators. FileNameJoin just assembles a file name; it does not search for the file specified.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (4)

```mathematica
In[1]:= FileNameJoin[{"usr", "local", "bin"}]
Out[1]= "usr/local/bin"
```

Components may themselves contain separators

```mathematica
In[2]:= FileNameJoin[{"a/b", "c"}]
Out[2]= "a/b/c"
```

An empty leading part gives an absolute path

```mathematica
In[3]:= FileNameJoin[{"", "etc", "hosts"}]
Out[3]= "/etc/hosts"
```

Backslash separator

```mathematica
In[4]:= FileNameJoin[{"dir1", "dir2"}, OperatingSystem -> "Windows"]
Out[4]= "dir1\dir2"
```

### Applications (4)

Join with the OS separator

```mathematica
In[5]:= FileNameJoin[{"usr", "local", "bin"}]
Out[5]= "usr/local/bin"
```

Components may themselves contain separators

```mathematica
In[6]:= FileNameJoin[{"a/b", "c"}]
Out[6]= "a/b/c"
```

An empty leading part gives an absolute path

```mathematica
In[7]:= FileNameJoin[{"", "etc", "hosts"}]
Out[7]= "/etc/hosts"
```

Backslash separator

```mathematica
In[8]:= FileNameJoin[{"dir1", "dir2"}, OperatingSystem -> "Windows"]
Out[8]= "dir1\dir2"
```

## Implementation notes

**Algorithm.** `builtin_filenamejoin` is a **pure string operation** — it never touches the
filesystem. It first decodes any trailing `OperatingSystem -> "Windows" | "MacOSX" | "Unix"`
option rules into a target separator (`sep`) and a `windows` flag (default: the host's
`HOST_SEP`, `/` on POSIX), rejecting an unknown OS or a non-option trailing argument by leaving
the call unevaluated. The components are then a single string (canonicalized) or a `List` of
strings; a non-string member returns `NULL`, and `FileNameJoin[]` prints `FileNameJoin::argx`.

`fnj_build` does the assembly. It emits a leading separator when the first component is empty or
separator-led (an **absolute path**), or two separators for a Windows UNC `\\server\share`
prefix. It then walks every component, skipping separator runs and copying each maximal
non-separator **segment**, inserting `sep` only between emitted segments. Empty segments from
leading, trailing, or duplicated separators are dropped, so `{"a//b", "c"}` collapses to
`"a/b/c"`.

**Data structures.** A small `const char**` of component pointers (a one-element stack buffer
for the lone-string case, else a `malloc`'d array borrowing the argument strings) and one
`malloc`'d output buffer sized from the summed component lengths. `fnj_is_sep` treats `/` as a
separator always and `\` as one only under `windows`.

**Complexity / limits.** `O(total input length)`, single pass. It is the inverse of
`FileNameSplit`: `FileNameJoin[FileNameSplit[name]]` reconstructs a
canonicalized `name`. `ATTR_PROTECTED`.

- `Protected`.
- Pure string operation — does not touch the filesystem.
- Components may themselves contain separators; each is split into segments and rejoined, so duplicate and trailing separators collapse (`{"a//b", "c"}` → `"a/b/c"`).
- An empty (or separator-led) leading component yields an absolute path: `{"", "usr", "bin"}` → `"/usr/bin"`.
- `"Windows"` uses `\` and preserves a leading `\\server\share` UNC prefix as a single unit; `"MacOSX"`/`"Unix"` use `/`. The default is the host operating system's separator.
- `Options[FileNameJoin]` reports the `OperatingSystem` default.
- `FileNameJoin[]` prints `FileNameJoin::argx` and stays unevaluated; a non-string/non-list argument, a list containing a non-string, or an unknown OS leaves the call unevaluated.

**Attributes:** `Protected`.

## References

- Source: [`src/files.c`](https://github.com/stblake/mathilda/blob/main/src/files.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
- Tests: [`tests/test_files.c`](https://github.com/stblake/mathilda/blob/main/tests/test_files.c)

## Notes & additional examples

### Notes

`FileNameJoin[{...}]` assembles a path from components, and `FileNameJoin["name"]`
canonicalizes a single name. It is a **pure string operation** — it never touches
the filesystem. Each component is split into segments and rejoined, so duplicate
and trailing separators collapse (`{"a//b", "c"}` becomes `"a/b/c"`), and an empty
(or separator-led) leading component yields an absolute path.

The separator defaults to the host operating system's; `OperatingSystem ->
"Windows" | "MacOSX" | "Unix"` selects it explicitly (`"Windows"` uses `\` and
preserves a leading UNC `\\server\share`). `FileNameJoin` inverts
`FileNameSplit`: `FileNameJoin[FileNameSplit[name]]`
reconstructs a canonicalized `name`.
