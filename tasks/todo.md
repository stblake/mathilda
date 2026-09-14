# Task: ReadList builtin in new src/io/ subsystem

## Plan
Filename-only ReadList with type-directed reading + separator options.

## Steps
- [x] 1. Scaffold + minimal `ReadList["file"]` — done (all types implemented together).
- [x] 2. All single-type readers + SYM_ names — done, smoke-tested at REPL.
- [x] 3. Group mode `{t1,...,tk}` (EndOfFile padding) + `n` limit — done.
- [x] 4. Options via options_extract + symtab_set_options — done, all 5 options verified.
- [x] 5. Tests: tests/test_readlist.c (27 cases) + CMake wiring — all pass.
- [x] 6. Docs (file-io.md, changelog, Mathilda_spec.md) — done. Audits: check-c99 PASS; valgrind = baseline (no ReadList leak). Packed/nd/compile audits running.

## Review

**Status: complete.** ReadList implemented in the new `src/io/` subsystem.

Files created: `src/io/readlist.{c,h}`, `tests/test_readlist.c`.
Files modified: `makefile` (SRC += io/*.c), `src/core.c` (readlist_init call),
`src/sym_names.{h,c}` (13 interned names), `src/info.c` (docstring),
`tests/CMakeLists.txt` (COMMON_SRC + readlist_tests target), `docs/spec/builtins/file-io.md`,
`docs/spec/changelog/2026-09-14.md`, `Mathilda_spec.md`.

Behaviour (filename-only; no stream layer):
- `ReadList["f"]` / `[..., type]` / `[..., {t1,..}]` (grouped, EndOfFile-padded) / `[..., types, n]`.
- Types: Byte, Character, Expression, Number, Real, Record, String, Word.
- Options via shared `options_extract` + `Options[ReadList]`: RecordSeparators,
  WordSeparators, TokenWords, NullRecords, NullWords (SetOptions works).
- Missing file → `ReadList::noopen` + $Failed; malformed number → `ReadList::readn` + $Failed slot.

Verification:
- Main binary builds clean (gcc-16, no warnings from readlist.c).
- `readlist_tests`: 27/27 pass. Neighbours (readwrite_tests, files_tests) still pass.
- valgrind: leak profile byte-identical to a no-ReadList baseline (no ReadList leak).
- Audits: check-c99 PASS; check-packed-aware / check-array-exactness / check-nd-surfaces PASS
  (ReadList correctly absent — not a numeric head).
- check-compile-coverage fails, but PRE-EXISTING on main (21 image/packed heads; identical list
  with/without my change, verified via a clean HEAD worktree). ReadList not in the list. Not a
  regression; left untouched (out of scope). See memory: project_check_compile_coverage_preexisting_red.

Not done (intentional): no OpenRead/Read/Close/Streams (no stream infra); no NDArray/Compile
fast paths (structural/I-O head, returns a List).
