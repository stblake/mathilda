---
source: src/files.c
---
`builtin_filebasename` is pure string manipulation (no filesystem access). It takes the filename component (everything after the last `/`, via `filename_component`), then strips a trailing extension using `extension_offset` — the suffix after the last `.` that is neither the first nor last character of the component. The base name is the component up to (not including) that `.`; with no qualifying extension the whole component (including any trailing `.`) is kept. The result is `memcpy`'d into a fresh buffer and returned as an `EXPR_STRING`. Non-string input is unevaluated. `ATTR_PROTECTED`.
