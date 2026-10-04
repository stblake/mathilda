---
source: src/files.c
---
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
