---
source: src/files.c
---
`builtin_fileextension` is pure string manipulation (no filesystem access). It isolates the filename component (after the last `/`, via `filename_component`) and finds the extension start with `extension_offset`: the offset just past the last `.` that is neither the first nor the last character of the component (so `.bashrc` and `file.` have no extension). It returns the substring from that offset — i.e. the suffix without the dot — as an `EXPR_STRING`, or `""` when there is no extension. `ATTR_PROTECTED`.
