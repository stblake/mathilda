# Fix passagemath aarch64-musl build break + warnings

CI: passagemath build on `aarch64-alpine-linux-musl`, GCC 14.2.0.
`make[3]: *** [src/print.o] Error 1`, `src/repl.o Error 1` — the build FAILS
(not just warns). Plus six -Wall/-Wextra warning sites.

## Root causes
- **ERROR (musl)**: musl declares `stdin/stdout/stderr` as `FILE *const`, so
  `stdout = …` / `stderr = …` is a hard error. glibc and macOS make them
  assignable, so this is invisible locally and to `make check-c99`.
  - `src/print.c` (expr_to_string / _fullform / numberform_format_result_to_string)
    swaps `stdout` to an `open_memstream` to capture its own output.
  - `src/repl.c` (pipe_capture_begin/end) swaps `stdout`+`stderr` to capture a
    script's output for the pipe protocol.
- **WARNING / latent bug (aarch64)**: `src/graph/galg_matching.c` uses
  `char* side` with sentinel `-1`; `char` is UNSIGNED on ARM, so `side[x] < 0`
  is always false and `>= 0` always true — bipartite detection is broken on
  aarch64 (flagged by -Wtype-limits).
- **WARNINGS (all platforms)**: unhandled enum in a switch, three
  misleading-indentation files, unused LAPACK-only locals.

## Plan
- [ ] print.h: expose a print-subsystem output sink (`mth_out`, push/pop).
- [ ] print.c: route printf/putchar/fputs through the sink; capture helpers
      push/pop the sink instead of assigning `stdout`.
- [ ] numberform.c: route its `fputs(…, stdout)` through `mth_out()` (same
      print subsystem, different TU).
- [ ] repl.c: rewrite pipe capture as fd-level redirect (tmpfile + dup2),
      portable on musl and a catch-all for all stdout/stderr.
- [ ] galg_matching.c: `char` → `signed char` for `side` (fixes warning AND the
      aarch64 correctness bug).
- [ ] integrate.c: handle `METHOD_INTEGRAL_REP` (definite-only → leave
      unevaluated) in the indefinite cascade switch.
- [ ] ndsolve_common.c: scope `n/nrhs/info` inside `#ifdef USE_LAPACK`.
- [ ] integrate_intrep.c / integrate_ramanujan.c / integrate_residue.c: split
      `if (p) free(p); return …;` one-liners so the `return` isn't mis-guarded.
- [ ] Add a `make check-c99` rule for stdin/stdout/stderr assignment (this class
      is invisible on glibc/macOS).
- [ ] version.h bump + changelog entry.
- [ ] Verify: local gcc build clean, `make check-c99`, targeted recompiles.

## Review

All items done; verified on GCC 16 (Homebrew), v0.342.

**Build break (the actual failure).**
- print.{c,h}: added a print-subsystem sink (`mth_out`/`mth_out_push`/
  `mth_out_pop`). print.c routes `printf`/`putchar` via file-local macros and
  `fputs(…, MTH_POUT)`; the three capture helpers push the sink to a memstream
  instead of assigning `stdout`. numberform.c's 15 `fputs(…, stdout)` now use
  `mth_out()`.
- repl.c + ffi/mathilda_ffi.c: pipe/notebook capture rewritten as fd-level
  redirect (`tmpfile` + `dup`/`dup2`), portable on musl and a catch-all. ffi had
  the *same* latent bug (CI died at print.o/repl.o before reaching it).
- No `stdin/stdout/stderr` assignment remains tree-wide (grep + new gate).

**Warnings.** integrate.c switch handles METHOD_INTEGRAL_REP; galg_matching.c
`char`→`signed char` (also fixes the aarch64 bipartite correctness bug);
ndsolve_common.c LAPACK locals scoped; three misleading-indentation files
reformatted.

**New gate.** tools/check_c99_portability.py now flags assignment to a standard
stream (invisible on glibc/macOS). `make check-c99` passes; self-tested to catch
a planted `stdout = m;` with no false positives (old_stdout=, ==stdout, ->stdout,
.stderr=).

**Verification.** Full tree builds warning-clean (USE_LAPACK on); ndsolve_common
clean without USE_LAPACK (the CI config); galg clean under `-funsigned-char`
(and old `char` confirmed to warn). Tests pass: print_tests, numberform_tests,
graph_algos_tests, graph_tests. `-file` printing and NDJSON pipe capture
(stdout stream + stderr message + multi-line) re-verified. check-messages OK.

**Not done (left to the user, per "commit only when asked").** No git commit.
When committing: message tag `; v0.342`, lightweight tag `v0.342` at that commit.
