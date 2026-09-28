# Task: Implement `DateList`

## Plan
- [x] 1. `src/datetime.c`: add `#include "message.h"`, `<stdarg.h>`; `dt_msg()` helper
- [x] 2. `src/datetime.c`: `gregorian_from_abstime()` inverse (JDN→Gregorian + intraday split)
- [x] 3. `src/datetime.c`: extract `datelist_parts_to_abstime()`, refactor `builtin_absolute_time` to use it
- [x] 4. `src/datetime.c`: string helpers — month/day-name tables, tokenizer, bare + format-form parse
- [x] 5. `src/datetime.c`: `builtin_date_list()` dispatcher (all forms)
- [x] 6. `src/datetime.c` `datetime_init()`: register + `ATTR_PROTECTED`
- [x] 7. `src/datetime.h`: prototype
- [x] 8. `src/sym_names.{h,c}`: `SYM_DateList` (decl + def + intern)
- [x] 9. `src/info.c`: docstring under `// Time and Date`
- [x] 10. `src/version.h`: bump 0.218 → 0.219 (number + string)
- [x] 11. `tests/test_datetime.c`: extend with DateList cases; `tests/CMakeLists.txt`: add `add_test(NAME datetime_tests ...)`
- [x] 12. `docs/spec/builtins/time-and-date.md` + `docs/spec/changelog/2026-09-28.md`
- [x] 13. Build main + tests; run all pasted examples; ctest; check-c99; leaks

## Decisions (from user)
- String parsing: **Full** (all pasted examples + common formats)
- TimeZone: **Defer** (local time, no TZ/DST correction — like AbsoluteTime)
- DateObject specs: out of scope (DateObject unimplemented)

## Review

Done. `DateList` implemented in `src/datetime.c`, v0.219.

**Design.** One backend: every form (now / absolute-time number / {y,m,...} spec /
"string" / {"string",{elems}}) builds a `parts[6]` array, converts to
seconds-since-1900 via the shared `datelist_parts_to_abstime` (factored out of
`builtin_absolute_time`, reusing `days_since_1900`), then inverts with the new
`gregorian_from_abstime` (inverse Fliegel & Van Flandern + intraday split). This
makes out-of-range/fractional fields normalize identically to AbsoluteTime and
Mathematica. Strings: `dt_tokenize` on non-alnum runs; bare parse handles month
names + 4-digit year + AM/PM, with US M/D/Y default + `DateList::ambig` for
order-ambiguous numeric strings and ISO for a leading 4-digit token; format form
ignores separator entries and assigns tokens to element entries in order.

**Verification.**
- All 18 pasted acceptance examples reproduce the expected Mathematica output
  (the {…8.1} residue lands on 8h 6m 0s vs MMA's 8h 6m 4.77e-7s — same minute, a
  harmless last-ULP difference).
- `datetime_tests`: all pass; now registered with ctest (`ctest -R datetime` →
  1/1 passed). It previously built but was never run.
- macOS `leaks --atExit`: **0 leaks / 0 bytes**.
- `make check-c99`: exit 0. Main build + relink clean under
  `gcc -std=c99 -Wall -Wextra` (+ the `-Werror=` set).
- `$VersionNumber` → 0.219; `Attributes[DateList]` → {Protected};
  `Information[DateList]` prints the docstring.

**Notes.**
- Removed a `NumericalOrder` mention from the docstring/spec: `NumericalOrder`
  is not a registered builtin in Mathilda (only referenced in an info.c
  docstring), so the claim would have been inaccurate.
- Packed/Compile surfaces: `DateList` is a structural date head taking a single
  spec (not an element-wise numeric kernel), so it is deliberately un-vectorized,
  exactly like its sibling `AbsoluteTime` (neither is on `AWARE`). It adds no
  NDArray dispatch site, so `check-packed-aware` and the curated-probe audits
  cannot newly flag it.
- Deferred (per user): `TimeZone`/`$TimeZone`; `DateObject` specs
  (`DateObject` unimplemented) are left unevaluated.
- Not yet committed/tagged — per project convention this substantive change
  should be committed with `; v0.219` and tagged `v0.219` when the user is ready.
