---
source: src/files.c
---
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
