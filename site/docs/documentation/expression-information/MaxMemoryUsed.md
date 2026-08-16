# MaxMemoryUsed

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MaxMemoryUsed[] gives the peak number of bytes resident for the Mathilda process over its lifetime. The high-water mark comes from the operating system, so it catches spikes that occurred between two calls to MemoryInUse rather than only the largest value previously observed.`**

## Examples

_No verified examples yet for this function._

## Algorithm

meminfo.c -- MemoryInUse and MaxMemoryUsed.

WHAT THESE MEASURE, AND HOW IT DIFFERS FROM MATHEMATICA. Wolfram's MemoryInUse[] reports the bytes used to store the data of the current session -- expressions, definitions, caches -- and nothing else. Mathilda reports the process's RESIDENT SET SIZE, which also includes the binary, the shared libraries (GMP, MPFR, LAPACK, Readline, and on a graphics build Raylib), the stacks, and whatever the allocator is holding but has not returned to the OS.

That is a real difference and it is documented rather than papered over, because the numbers are not interchangeable: on a freshly started kernel Wolfram's figure is small while this one is tens of megabytes of mapped libraries. Reporting session-data bytes exactly would mean routing every allocation in the tree through a counting wrapper -- about 500 modules, all of them calling malloc directly today -- and the resulting number would still miss what the allocator retains.

RSS is also the more USEFUL quantity for the thing that prompted this: a notebook status bar. It is the number Activity Monitor and top show, so a user comparing the two sees them agree, where a smaller session-only figure would look like an under-report.

PORTABILITY, which has two traps in it.

```text
`getrusage` is POSIX, not C99, so the feature-test macro has to come BEFORE any include --
```

below the first one the header has already been parsed with the wrong namespace and the macro does nothing. glibc hides the symbol under -std=c99 while Darwin exposes it anyway, so getting this wrong compiles clean here and breaks only on Linux. Same pattern as src/ndkernels.c, src/core.c, src/repl.c.

```text
`ru_maxrss` IS NOT IN THE SAME UNIT ON BOTH PLATFORMS. Darwin reports BYTES; Linux
```

reports KILOBYTES. A single unconverted use is wrong by a factor of 1024 on one of the two and silently plausible on both -- exactly the class of bug the int64_t/long long note in CLAUDE.md is about. Hence the explicit per-platform scaling below.

AND THE FEATURE-TEST MACRO RUNS THE OPPOSITE WAY FROM EVERY OTHER FILE IN THE TREE, which is why the guard below is conditional rather than the usual unconditional _XOPEN_SOURCE.

```text
`ru_maxrss` is a BSD extension, not POSIX -- POSIX only requires ru_utime and ru_stime. So
```

on Darwin, defining _XOPEN_SOURCE selects the strict POSIX subset and HIDES the field:

```text
`no member named 'ru_maxrss' in 'struct rusage'`, from the very macro that every other
```

module needs in order to see POSIX at all. Meanwhile glibc needs a feature macro to expose

```text
`getrusage` itself. Each platform therefore gets the macro that widens ITS namespace --
```

_DARWIN_C_SOURCE for the BSD fields, _DEFAULT_SOURCE plus _XOPEN_SOURCE for glibc -- and both still come before any include, where a feature-test macro has to be to do anything.

## Implementation notes

- `Protected`.
- A genuine high-water mark from the operating system (`getrusage`'s `ru_maxrss`), not the
  largest value some earlier call to `MemoryInUse` happened to observe. That distinction is
  the point: a polled maximum misses any spike falling between two polls, and a status bar
  polling once a second would miss nearly every spike worth knowing about.
- `MaxMemoryUsed[] >= MemoryInUse[]` always holds, and the peak never decreases.
- **`ru_maxrss` is in different units on the two platforms** — bytes on Darwin, kilobytes on
  Linux — so it is scaled per platform. An unconverted use is wrong by a factor of 1024 on
  one of them while looking plausible on both.
- Its feature-test guard runs **opposite** to every other file in the tree. `ru_maxrss` is a
  BSD extension rather than POSIX (POSIX requires only `ru_utime` and `ru_stime`), so on
  Darwin defining `_XOPEN_SOURCE` selects the strict POSIX subset and *hides* the field —
  from the very macro every other module needs in order to see POSIX at all. So each platform
  gets the macro that widens *its* namespace: `_DARWIN_C_SOURCE` on macOS, `_DEFAULT_SOURCE`
  plus `_XOPEN_SOURCE` on glibc.

**Attributes:** `Protected`.

## References

**See also:** [MemoryInUse](../../expression-information/MemoryInUse/)

- Source: [`src/meminfo.c`](https://github.com/stblake/mathilda/blob/main/src/meminfo.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_meminfo.c`](https://github.com/stblake/mathilda/blob/main/tests/test_meminfo.c)
