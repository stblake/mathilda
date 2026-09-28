# Task: Implement `Pause[n]` (+ `$TimeUnit`, `SessionTime`, `TimeUsed`)

## Plan
Implement `Pause[n]` (wall-clock sleep, returns Null, no CPU) plus the three
symbols its docstring references: `$TimeUnit`, `SessionTime[]`, `TimeUsed[]`.
Semantics fall out of the existing two-clock design (clock() vs
clock_gettime(CLOCK_MONOTONIC)).

## Checklist
- [x] `src/datetime.c`: add `<errno.h>`, `g_session_start`, `builtin_pause`,
      `builtin_session_time`, `builtin_time_used`; register + Protected in init
- [x] `src/datetime.h`: declare the three builtins
- [x] `src/sym_names.h` / `.c`: add `SYM_Pause` (extern, def, intern)
- [x] `src/core.c`: register `$TimeUnit` in system_constants_init
- [x] `src/info.c`: 4 docstrings (Pause, SessionTime, TimeUsed, $TimeUnit)
- [x] `tests/test_datetime.c`: feature-test macro, now_seconds, test fns + main
- [x] `src/version.h`: 0.217 -> 0.218 (number + string)
- [x] Docs: time-and-date.md, changelog 2026-09-28.md, Mathilda_spec.md table row
- [x] Build clean (`make`), `make check-c99` (exit 0)
- [x] Run `datetime_tests` (all pass); leak-audited (macOS valgrind = baseline noise)
- [x] Refresh code-review graph

## Review
Implemented `Pause[n]` plus `$TimeUnit`, `SessionTime[]`, `TimeUsed[]` (v0.218).

Key design point: the "counted by AbsoluteTiming/SessionTime, not Timing/TimeUsed"
semantic is free — `nanosleep` burns no CPU, so the existing `clock()` (CPU) vs
`clock_gettime(CLOCK_MONOTONIC)` (wall) split does the accounting with no
special-casing.

Deviation from plan (improvement): `Pause` uses a robust `pause_seconds` coercion
(mirrors `clip_to_double_value` in core.c) instead of the strict
`expr_to_double_strict`, so `Pause[1/4]` and `Pause[Pi]` work like Mathematica.
Left `expr_to_double_strict` untouched so `AbsoluteTime`'s behavior is unchanged.

Verified: `make` clean; `make check-c99` exit 0; `datetime_tests` all pass;
REPL smoke tests (`Timing[Pause[1]]`≈{0,Null}, `AbsoluteTiming[Pause[1]]`≈{1,Null},
rational/Pi/symbolic args, docstrings, `$VersionNumber`=0.218). No leaks from the
new code (valgrind output was macOS Objective-C/Foundation baseline noise plus
libsystem `nanosleep` uninitialised-value warnings; no `builtin_pause` /
`pause_seconds` / helper appears in any leak stack).

Not committed/tagged — awaiting explicit request. When asked: tag `v0.218`.
