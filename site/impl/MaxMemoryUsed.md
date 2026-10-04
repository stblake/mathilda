---
source: src/meminfo.c
---
**Algorithm.** `builtin_maxmemoryused` (`src/meminfo.c`) returns the peak resident
bytes over the life of the process, read by `meminfo_peak` from `getrusage`'s
`ru_maxrss`. It is a genuine OS high-water mark, **not** the largest value some
earlier `MemoryInUse[]` call happened to observe: a polled maximum misses any
spike falling between two polls, and a once-a-second status bar would miss nearly
every spike worth knowing about. The peak and the current RSS come from two
different OS counters (`getrusage` here versus the per-task resident size in
`MemoryInUse[]`), so the two are not guaranteed to agree to the byte at a given
instant — near startup the `ru_maxrss` peak can read a page below the mach current
RSS.

**Data structures.** None of its own — one `getrusage(RUSAGE_SELF, ...)` call
filling a `uint64_t`, returned as an `EXPR_INTEGER`. Takes no arguments; a
non-empty call returns `NULL` (stays unevaluated).

**Complexity / limits.** `O(1)`. Two portability traps are handled explicitly.
`ru_maxrss` is **bytes on Darwin but kilobytes on Linux**, so the Linux path
scales by 1024 — an unconverted use is wrong by that factor on one platform while
looking plausible on both. And its feature-test guard runs *opposite* to every
other file: `ru_maxrss` is a BSD extension (not POSIX), so `_XOPEN_SOURCE` would
*hide* the field on Darwin; the file defines `_DARWIN_C_SOURCE` on macOS and
`_DEFAULT_SOURCE` + `_XOPEN_SOURCE` on glibc, each before any `#include`. **The
value changes from run to run**, so examples assert only structure
(`Head[MaxMemoryUsed[]]` is `Integer`, `MaxMemoryUsed[] > 0`), never a literal byte
count. Attributes `Protected`.
