---
source: src/datetime.c
---
**Algorithm.** `builtin_date_list` reduces every surface form — the current time
(`DateList[]`), an absolute-time number (seconds since 1900-01-01), a
`{y,m,d,h,mi,s}` spec with right-elision, a free-form date string, or a
`{"string", {elements}}` format spec — to a raw `parts[6]` array, then runs it
through one shared pipeline: `datelist_parts_to_abstime` →
`gregorian_from_abstime` → `datelist_make_result`. So out-of-range fields
normalise identically to `AbsoluteTime` and to Mathematica
(`DateList[{2022,2,31}]` is `{2022,3,3,0,0,0.}`), and `DateList[AbsoluteTime[s]]`
round-trips. The result is five exact integers plus a `Real` second, so a
current-time or fractional result carries sub-second precision.

**Calendar arithmetic.** `datelist_parts_to_abstime` splits the day into an
integer part (fed to `days_since_1900`, the Fliegel & Van Flandern Julian-Day
formula, offset to the 1900 epoch) and a fractional part carried as seconds;
year and month must be integer-valued (the lengths of years and months vary) and
a non-integer month gives `DateList::arg`. `gregorian_from_abstime` is the exact
integer inverse. The string path (`datelist_from_string`) tokenises on any
non-alphanumeric run, recognises month/weekday names and am/pm, and warns
`DateList::ambig` on an order-ambiguous all-numeric string (US `M/D/Y` default);
the format path (`datelist_from_format`) assigns tokens to element codes in order.

**Data structures / limits.** `double parts[6]` and fixed `char[16][64]` token
buffers; all calendar math is `int64_t`. Local time, with no time-zone / DST /
leap-second correction. `ATTR_PROTECTED`. Structural calendar logic, not a
numeric buffer op — no packed/`Compile[]` surface.
