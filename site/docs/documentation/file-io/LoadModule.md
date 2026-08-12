# LoadModule

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`LoadModule["relpath"]`**

loads the internal Mathilda source module at relpath (relative to src/internal), resolving the location independently of the current working directory. Each module is loaded at most once. Returns True if the module was located and loaded, False otherwise.

## Examples

_No verified examples yet for this function._

## Algorithm

loadmodule.c -- runtime loading of Mathilda (.m) source modules.

See loadmodule.h for the contract. The file-reading loop here is the same one Get[] uses (builtin_get in readwrite.c is now a thin wrapper over mathilda_run_file); mathilda_load_module adds working-directory-independent path resolution and load-once bookkeeping on top.

Path resolution (mathilda_resolve_internal) is deliberately CWD-independent: a relocated or installed binary must still find its bundled src/internal tree. On glibc, readlink("/proc/self/exe") needs the POSIX feature-test macro to be visible under -std=c99; define it before any system header.

## Implementation notes

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

## See also

[Get](../../file-io/Get/)

## References

- Source: [`src/loadmodule.c`](https://github.com/stblake/mathilda/blob/main/src/loadmodule.c)
- Specification: [`docs/spec/builtins/file-io.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/file-io.md)
