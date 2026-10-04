---
source: src/loadmodule.c
---
**Algorithm.** `builtin_loadmodule` takes one `EXPR_STRING` relpath and calls
`mathilda_load_module`, returning `True`/`False`. `mathilda_load_module` short-circuits to
`True` if the relpath is already in the load-once table (`lm_already_loaded`); otherwise it
resolves the path, runs the file, and on success records the relpath (`lm_mark_loaded`) so a
repeat call never re-reads it. A relpath that resolves nowhere emits `LoadModule::nofile`
through `mth_message` (the `Quiet`/`Check` funnel) and returns `False`.

`mathilda_resolve_internal` is deliberately **working-directory-independent**, trying in order:
(1) `$MATHILDA_HOME/<relpath>`; (2–3) relative to the running executable —
`<exe_dir>/src/internal/<relpath>` and `<exe_dir>/../share/mathilda/internal/<relpath>`, via
`readlink("/proc/self/exe")` on Linux, `_NSGetExecutablePath` on macOS, `GetModuleFileNameA` on
Windows — so a relocated or installed binary still finds its bundled modules; (4) a compile-time
`MATHILDA_PREFIX/share/mathilda/internal/`; (5) a CWD ladder (`src/internal/` up to three levels
up). The winning base directory is cached (`lm_base`) so later lookups skip the search.

`mathilda_run_file` is the file-reading core `Get` now wraps: slurp the file into a
buffer, then loop `parse_next_expression(&ptr)` / `evaluate` over it, freeing each parse and
keeping the last evaluated value (`Null` for an empty file). `LoadModule` discards that value
and reports only whether the file was opened.

**Data structures.** Fixed-size bookkeeping tables (`lm_loaded[256][256]` relpaths; a cached
`lm_base`) and a transient slurp buffer. On glibc `_POSIX_C_SOURCE 200112L` is defined before
any include so `readlink` is visible under `-std=c99` (SPEC.md §10).

**Complexity / limits.** First load is `O(file size)` to read and evaluate; a repeat is `O(1)`
(table hit). Each module is loaded **at most once**, so the lazy per-family loading in
`FullSimplify` never re-registers its rules. `ATTR_PROTECTED`.
