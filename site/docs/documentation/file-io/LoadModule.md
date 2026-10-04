# LoadModule

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LoadModule["relpath"]`**

loads the internal Mathilda source module at relpath (relative to src/internal), resolving the location independently of the current working directory. Each module is loaded at most once. Returns True if the module was located and loaded, False otherwise.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic Examples (3)

Load an internal module

```mathematica
In[1]:= LoadModule["simp/FullSimplify.m"]
Out[1]= True
```

Already loaded: still True, not re-read

```mathematica
In[2]:= LoadModule["simp/FullSimplify.m"]
Out[2]= True
```

Not found -> False

```mathematica
In[3]:= Quiet[LoadModule["no/such/module.m"]]
Out[3]= False
```

### Applications (3)

Load an internal module; True on success

```mathematica
In[4]:= LoadModule["simp/FullSimplify.m"]
Out[4]= True
```

Already loaded: still True, not re-read

```mathematica
In[5]:= LoadModule["simp/FullSimplify.m"]; LoadModule["simp/FullSimplify.m"]
Out[5]= True
```

Not found -> False

```mathematica
In[6]:= Quiet[LoadModule["no/such/module.m"]]
Out[6]= False
```

## Algorithm

loadmodule.c -- runtime loading of Mathilda (.m) source modules.

See loadmodule.h for the contract. The file-reading loop here is the same one Get[] uses (builtin_get in readwrite.c is now a thin wrapper over mathilda_run_file); mathilda_load_module adds working-directory-independent path resolution and load-once bookkeeping on top.

Path resolution (mathilda_resolve_internal) is deliberately CWD-independent: a relocated or installed binary must still find its bundled src/internal tree. On glibc, readlink("/proc/self/exe") needs the POSIX feature-test macro to be visible under -std=c99; define it before any system header.

## Implementation notes

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

- `Protected`. Returns `True` if the module was located and loaded (or had already
  been loaded), `False` otherwise.
- Resolution is independent of the current working directory and tries, in
  order: `$MATHILDA_HOME/<relpath>`; `<exe_dir>/src/internal/<relpath>` and
  `<exe_dir>/../share/mathilda/internal/<relpath>` (relative to the running
  binary, so a relocated or installed executable still finds its modules);
  `$(PREFIX)/share/mathilda/internal/<relpath>` when built with a compile-time
  `MATHILDA_PREFIX`; and finally a CWD ladder (`src/internal/`,
  `../src/internal/`, `../../src/internal/`, `../../../src/internal/`). This
  works from the REPL (run at the repo root or from anywhere with the binary in
  place), from the test binaries (run from `tests/build/`), and from a binary
  copied to a `bin` directory with `MATHILDA_HOME` pointing at `src/internal`.
  The winning base directory is cached after the first successful lookup.
- Each module is loaded **at most once**, so repeated calls — and the lazy
  per-family loading used by [`FullSimplify`](../simplification/index.md) —
  never re-register rules.
- Generalises the bespoke fallback previously hard-coded for the CRC integral
  tables; `Get` (above) shares its file-reading core.

**Attributes:** `Protected`.

## References

**See also:** [Get](../../file-io/Get/)

- Source: [`src/loadmodule.c`](https://github.com/stblake/mathilda/blob/main/src/loadmodule.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)

## Notes & additional examples

### Notes

`LoadModule["relpath"]` loads a Mathilda `.m` source module, with `relpath`
relative to the source tree's `src/internal` directory (e.g.
`"simp/FullSimplify.m"`). It returns `True` if the module was found and loaded (or
had already been loaded) and `False` otherwise. This is the mechanism the lazy
per-family loading in `FullSimplify` uses.

Resolution is deliberately **independent of the working directory**: it tries
`$MATHILDA_HOME`, then paths relative to the running executable (so a relocated or
installed binary still finds its bundled modules), then a compile-time prefix, then
a CWD ladder. Each module is loaded **at most once**, so repeated calls never
re-register rules — which is why the second call above is still `True` without
re-reading the file. It shares its file-reading core with `Get`.
