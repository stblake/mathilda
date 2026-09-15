# Task: Read + a first-class stream layer for src/io/

## Plan
Introduce Mathilda's first input/output stream layer. `Read` is the single-object
primitive with a persistent current point; `ReadList` is refactored to loop the
same shared reading engine. Input + output + inspection/positioning.

Builtins: Read, OpenRead, OpenWrite, OpenAppend, Write, WriteString, Close,
Streams, StreamPosition, SetStreamPosition. Inert heads: InputStream,
OutputStream, File.

## Steps
- [x] 1. Stream registry + objects + OpenRead/OpenWrite/OpenAppend/Close/Streams — done.
- [x] 2. Reading engine (read.c) + read_structure (nested + EndOfFile pad, unevaluated Expression leaf) + builtin_read + Options[Read] — done.
- [x] 3. Refactor ReadList onto the shared engine — done, readlist_tests 27/27 green.
- [x] 4. Write/WriteString + StreamPosition/SetStreamPosition — done.
- [x] 5. tests/test_streams.c (36 cases) + CMake wiring — all pass.
- [x] 6. Docs + audits — done. check-c99 PASS; valgrind clean (no Mathilda frame in any leak).

## Review

**Status: complete.** First input/output stream layer added to `src/io/`.

Files created: `src/io/read.{c,h}` (shared reading engine + Read), `src/io/streams.{c,h}`
(registry + management/output builtins + inert heads), `tests/test_streams.c`.
Files modified: `src/io/readlist.c` (refactored onto the shared engine),
`src/sym_names.{h,c}` (13 interned names), `src/core.c` (streams_init/read_init),
`src/info.c` (10 docstrings), `tests/CMakeLists.txt` (COMMON_SRC + streams_tests),
`docs/spec/builtins/file-io.md`, `docs/spec/changelog/2026-09-14.md`, `Mathilda_spec.md`.

Behaviour:
- `Read[stream|"file"|File["file"], spec?]` — one object; nested type structures
  (depth-first), `Hold[Expression]` held vs `Expression` evaluated, arbitrary heads.
  EndOfFile past EOF; $Failed for closed stream / malformed number.
- Persistent current point: successive Read/ReadList advance; filename form auto-opens
  and stays open; `OpenRead`/`OpenWrite`/`OpenAppend` → InputStream/OutputStream handles.
- `Write`/`WriteString` (fflush, round-trips with Read/ReadList); `Close` (returns name);
  `Streams[]`/`Streams["file"]`; `StreamPosition`/`SetStreamPosition[..., n|Infinity]`.
- `ReadList` is now literally a loop of the shared single-object reader; also reads from
  an open InputStream.

Verification:
- Main binary builds clean (gcc-16, no warnings from read.c/streams.c/readlist.c).
- streams_tests: 36/36 pass. readlist_tests: 27/27 pass (behaviour-preserving refactor).
- valgrind on streams_tests: 0 leak contexts touch any Mathilda frame (all "definitely
  lost" blocks are macOS objc/dyld baseline noise). Registry freed by Close + atexit.
- Audits: check-c99 PASS; check-packed-aware PASS (782 builtins, Read/stream heads
  correctly NOT on any NDArray surface — they return objects/lists, not numeric kernels).
- check-compile-coverage remains the PRE-EXISTING image-head failure (Read/stream heads
  absent from it); not a regression. See memory project_check_compile_coverage_preexisting_red.

Not done (intentional/out of scope): ReadString/ReadLine/Skip/Find; pipes; OpenRead
options stored on the stream; $Input/$Output redirection; exact syntax-error point (a
parse failure stops, as ReadList already does).
