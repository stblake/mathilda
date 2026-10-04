---
source: src/datetime.c
---
**Algorithm.** `builtin_unix_time` gives seconds since the Unix epoch
1970-01-01 00:00:00 GMT. `UnixTime[]` returns `time(NULL)` directly — the true
POSIX second. Every dated form (a number taken as an `AbsoluteTime`, a
`{y,m,d,h,mi,s}` spec with elision, a date string, or a `{"string",{elements}}`
spec) reuses the `DateList` backend to reach absolute seconds-since-1900, then
`unixtime_from_abstime` subtracts the fixed offset
`DT_UNIX_EPOCH_OFFSET = 2208988800` (= `days_since_1900(1970,1,1) * 86400`) and
rounds to the nearest whole second. Hence `UnixTime[spec] == AbsoluteTime[spec]
- 2208988800`, `UnixTime[{1970,1,1,0,0,0}] == 0`, and the string/format forms
interpret a spec exactly as `DateList` does.

**Integer result.** Unlike `AbsoluteTime`, `UnixTime` always returns an
`EXPR_INTEGER` — a fractional-second spec is rounded (`floor(x + 0.5)`) rather
than kept as a `Real`, a deliberate match to Mathematica. A `Real` comes back
only in the corner case where the value does not fit an `int64`. A non-integer
year or month gives `UnixTime::arg` and leaves the call unevaluated.

**Data structures / limits.** Shares `datelist_from_string` /
`datelist_from_format` / `datelist_result_to_abstime` with `DateList`, so the two
heads never disagree on how a spec is read. `time_t`/`int64_t`; local calendar
arithmetic for dated forms, GMT for `UnixTime[]`; no TZ/DST/leap-second
correction. `ATTR_PROTECTED`. Scalar calendar arithmetic, no packed/`Compile[]`
surface.
