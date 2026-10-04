---
source: src/meminfo.c
---
**Algorithm.** `builtin_memoryinuse` (`src/meminfo.c`) returns the process's
current **resident set size** in bytes, read by `meminfo_current`:
`mach_task_basic_info`'s `resident_size` on macOS (already bytes), and the
resident-pages field of `/proc/self/statm` scaled by `sysconf(_SC_PAGESIZE)` on
Linux (not an assumed 4096, so a 16 KiB-page kernel is correct). This is
deliberately **not** Mathematica's quantity — Wolfram counts only the current
session's data, whereas RSS also includes the binary, the shared libraries (GMP,
MPFR, LAPACK, Readline, Raylib), the stacks, and allocator slack. RSS was chosen
because it is the number Activity Monitor and `top` report, which is what a
notebook status bar wants to agree with.

**Data structures.** None of its own — one OS query filling a `uint64_t`, returned
as an `EXPR_INTEGER`. `MemoryInUse[]` takes no arguments; the one-argument
subkernel form is rejected (`arg_count != 0` returns `NULL`) because there are no
subkernels to ignore silently.

**Complexity / limits.** `O(1)` (a syscall / small file read). On a platform that
offers no way to ask, `meminfo_current` fails and the builtin returns `NULL`
(stays unevaluated) rather than reporting `0` — a zero would read as "no memory
in use", which is false and plausible in a status bar. **The value changes from
run to run**, so examples assert only structure (`Head[MemoryInUse[]]` is
`Integer`, `MemoryInUse[] > 0`), never a literal byte count. Attributes
`Protected`.
