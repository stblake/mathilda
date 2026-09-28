# Task: Implement `UnixTime`

## Plan
- [x] 1. `src/datetime.c`: `datelist_result_to_abstime()` helper (list → seconds-since-1900)
- [x] 2. `src/datetime.c`: `builtin_unix_time()` dispatcher (now / number / {y,m,...} / "string" / {"string",{elems}})
- [x] 3. `src/datetime.c` `datetime_init()`: register + `ATTR_PROTECTED`
- [x] 4. `src/datetime.h`: prototype
- [x] 5. `src/sym_names.{h,c}`: `SYM_UnixTime` (decl + def + intern)
- [x] 6. `src/info.c`: docstring under `// Time and Date`
- [x] 7. `src/version.h`: bump 0.219 → 0.220 (number + string)
- [x] 8. `tests/test_datetime.c`: extend with UnixTime cases (no CMake change needed)
- [x] 9. `docs/spec/builtins/time-and-date.md` + `docs/spec/changelog/2026-09-28.md`
- [x] 10. Build main + tests; run datetime_tests; REPL smoke; check-c99; leaks

## Decisions (from user + module conventions)
- `UnixTime[]`: **true POSIX `time(NULL)`** (real GMT epoch second, matches `date +%s`)
- Epoch offset: `days_since_1900(1970,1,1) * 86400 = 2208988800`
- Return type: **always an Integer** (nearest whole second, `floor(unix+0.5)`),
  MMA-faithful; Real only as an int64 overflow fallback. Differs from AbsoluteTime.
- Numeric arg = AbsoluteTime spec (seconds since 1900): `UnixTime[t] = t - offset`
- String forms reuse DateList's parser verbatim (via new list→abstime helper)
- No NDArray/packed/Compile surfaces (scalar date head, like AbsoluteTime/DateList)

## Review

Done. `UnixTime` implemented in `src/datetime.c`, v0.220.

**Design.** `UnixTime` is `AbsoluteTime` shifted by the fixed 1900→1970 epoch
offset (`2208988800`), always rounded to a whole-second integer. `UnixTime[]`
uses `time(NULL)` (true POSIX GMT second, matching `date +%s`); the numeric /
date-list / string forms reuse the DateList backend so the two heads agree on
interpretation. One new helper `datelist_result_to_abstime` reads a DateList
result list back to seconds-since-1900, letting the string/format forms reuse
the existing (tested) parser verbatim with no refactor.

**Verification.**
- `make` — main binary links cleanly with `datetime.o`; banner reports 0.220.
- REPL smoke (`-file`): epoch-zero=0, 2022=1640995200, elision, numeric arg,
  AbsoluteTime relation=True, string/format parity, fractional rounding
  (.4→…200, .6→…201), normalization, unevaluated `UnixTime[x]`, `{Protected}`.
- `datetime_tests` — all pass (12 new UnixTime cases).
- `make check-c99` — PASS (no new POSIX symbols; `time()` is C89).
- valgrind — no leak/error stack touches `builtin_unix_time`,
  `unixtime_from_abstime`, or `datelist_result_to_abstime`. Datetime-adjacent
  valgrind chatter is pre-existing macOS `libsystem` noise (`localtime` tz init
  in `AbsoluteTime`, `nanosleep` in `Pause`) plus ObjC/dispatch baseline —
  documented macOS valgrind behaviour, not real leaks.

**Not committed** — left for the user. Suggested: commit ending `; v0.220` and
tag `v0.220` (`git tag v0.220 && git push --follow-tags`).
