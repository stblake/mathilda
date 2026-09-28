# Task: Implement `FileSize`

## Plan
- [x] 1. `src/files.c`: add `#include "message.h"` + `<stdarg.h>`
- [x] 2. `src/files.c`: static `fs_msg()` Quiet/Check-aware message helper (mirror `dt_msg`)
- [x] 3. `src/files.c`: `builtin_filesize()` — stat → st_size Integer; missing → msg + $Failed
- [x] 4. `src/files.c` `files_init()`: register `FileSize` + `ATTR_PROTECTED`
- [x] 5. `src/files.c` / `src/files.h`: update header comment; add prototype
- [x] 6. `src/sym_names.{h,c}`: `SYM_FileSize` (decl + def + intern)
- [x] 7. `src/info.c`: docstring after `FilePrint` (no examples)
- [x] 8. `src/version.h`: bump 0.220 → 0.221 (number + string)
- [x] 9. `docs/spec/builtins/file-io.md` (inventory + section) + `docs/spec/changelog/2026-09-28.md`
- [x] 10. Build main; check-c99; REPL smoke; `tests/test_files.c` (7 cases); rebuild graph

## Decisions (from user + module conventions)
- Returns a plain **Integer** byte count — explicitly NOT a Quantity
- Missing file → `FileSize::nffil` via **Quiet/Check-aware** helper (user choice), return `$Failed`
- Uses `stat()` (follows symlinks); non-string / wrong arity → `NULL` (unevaluated), like `FileExistsQ`
- No NDArray/packed/Compile surfaces (structural IO head, scalar result)

## Review

Done. `FileSize` implemented in `src/files.c`, v0.221.

**Design.** `FileSize["name"]` mirrors its sibling `FileExistsQ`: one string arg,
`stat()` (follows symlinks), returns `expr_new_integer((int64_t)st.st_size)` — a
plain Integer, deliberately not a Quantity. A path that cannot be stat'd emits
`FileSize::nffil` through a new static `fs_msg` helper (copied from datetime.c's
`dt_msg`: `mth_msg_note_fired()` then early-return on `mth_msg_suppressed()`) and
returns `$Failed`, so the message respects `Quiet[]` and is visible to `Check[]`.
Wrong arity / non-string arg → `NULL` (unevaluated), so symbolic args flow
through. Not a numeric array kernel, so no packed/NDArray/`Compile[]` surfaces
(consistent with the other `File*` heads).

**Verification.**
- `make` — links cleanly with GCC 16; banner reports 0.221.
- REPL smoke (`-file`): `FileSize["src/files.c"]` = 32130 = `wc -c` exactly;
  `Head` → `Integer`; `Attributes` → `{Protected}`; missing → message + `$Failed`;
  `Quiet[...]` → `$Failed` with no message; `Check[...]` → caught; symbolic and
  wrong-arity left unevaluated.
- `make check-c99` — PASS (stat already guarded by `_POSIX_C_SOURCE`; no new
  POSIX symbol).
- `tests/test_files.c` — 7 new `FileSize` cases (exact bytes, empty→0, Integer
  head, missing→`$Failed` under Quiet, Check detection, bad-args, Protected); all
  85 tests pass, no `FileSize::nffil` leak to stderr (confirms Quiet path).
- Code-review graph rebuilt after edits.

**Not committed** — left for the user. Suggested: commit ending `; v0.221` and
tag `v0.221` (`git tag v0.221 && git push --follow-tags`).
