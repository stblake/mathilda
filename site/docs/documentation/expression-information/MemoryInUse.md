# MemoryInUse

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MemoryInUse[] gives the number of bytes of memory currently resident for the Mathilda process. This is the process resident set size, so unlike Mathematica's MemoryInUse it also counts the binary, the shared libraries and whatever the allocator holds without returning it to the system -- it is the figure Activity Monitor and top report, not a count of session data alone. Returns unevaluated on a platform that offers no way to ask, rather than reporting zero.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= MemoryInUse[]
Out[1]= 7249920

In[2]:= N[MemoryInUse[]/1024^2, 4]
Out[2]= 6.9297
```

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
- **This is the process resident set size, which is not the same quantity Mathematica
  reports.** Wolfram's `MemoryInUse[]` counts the bytes holding the current session's data —
  expressions, definitions, caches — and nothing else. Mathilda's also includes the binary,
  the shared libraries (GMP, MPFR, LAPACK, Readline, and Raylib on a graphics build), the
  stacks, and whatever the allocator is holding without having returned it to the system. On
  a freshly started kernel Wolfram's figure is small where this one is tens of megabytes of
  mapped libraries, so the two are **not interchangeable**.
- The difference is deliberate on both counts. Reporting session-data bytes exactly would
  mean routing every allocation in the tree through a counting wrapper — around 500 modules,
  all calling `malloc` directly today — and the result would still miss what the allocator
  retains. And RSS is the more useful number for the purpose that prompted this, a notebook
  status bar: it is what Activity Monitor and `top` show, so a user comparing them sees them
  agree.
- Reads `mach_task_basic_info` on macOS and `/proc/self/statm` on Linux, scaled by the actual
  page size rather than an assumed 4096.
- **Returns unevaluated** on a platform offering no way to ask, rather than reporting `0` — a
  zero would read as "no memory in use", which is false and looks entirely plausible in a
  status bar.
- `MemoryInUse[subkernel]` is not supported: there are no subkernels, so an argument returns
  unevaluated rather than being accepted and ignored.

**Attributes:** `Protected`.

## References

- Source: [`src/meminfo.c`](https://github.com/stblake/mathilda/blob/main/src/meminfo.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_meminfo.c`](https://github.com/stblake/mathilda/blob/main/tests/test_meminfo.c)
