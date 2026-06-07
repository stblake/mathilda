# AbsoluteTime

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AbsoluteTime[]
    gives the total number of seconds since the beginning of January 1, 1900, in your time zone.
AbsoluteTime[date]
    gives the absolute time specification corresponding to the given date specification.

The supported date specifications are:
    {y, m, d, h, m, s}    DateList specification
    time            AbsoluteTime specification (a number, returned unchanged)

DateList entries may be elided from the right: {y}, {y, m}, {y, m, d}, etc. fill the
missing fields with {_, 1, 1, 0, 0, 0}. Day, hour, minute, and second values may be
noninteger; the year and month must be integers. Date lists are converted to standard
normalized form, so e.g. AbsoluteTime[{2022, 2, 31}] = AbsoluteTime[{2022, 3, 3}].

AbsoluteTime[] uses whatever date and time have been set on your computer system. It
performs no corrections for time zones, daylight saving time, or leap seconds.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_absolute_time` returns seconds since the Mathematica epoch 1900-01-01. With no argument it reads `time(NULL)` + `localtime`, splits the broken-down time into y/m/d/h/m/s, and combines `days_since_1900(...)*86400 + h*3600 + m*60 + s` (returned as `EXPR_REAL`). A numeric argument passes through as an existing absolute-time spec. A `{y,m,d,h,m,s}` `List` (1–6 elements, defaults `{y,1,1,0,0,0}`) is converted the same way; year and month must be integral, while day/hour/minute/second may be fractional (the day's fractional part is folded into seconds). When all inputs are integral and the total is representable, it returns an exact `EXPR_INTEGER`, otherwise an `EXPR_REAL`.

**Data structures.** Day numbers come from `days_since_1900`, which applies the Fliegel–Van Flandern Julian-Day-Number formula (offset by the JDN of 1900-01-01) and normalises out-of-range months explicitly. `ATTR_PROTECTED`.

- `Protected`.
- Year and month must be integer-valued; day, hour, minute, and second may be noninteger.
- Out-of-range date components are converted to standard normalized form, e.g. `AbsoluteTime[{2022, 2, 31}] == AbsoluteTime[{2022, 3, 3}] == 3855254400`.
- Performs no corrections for time zones, daylight saving time, or leap seconds.
- Returns an integer when every component is integer-valued and the total is exact; otherwise returns a real.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/datetime.c`](https://github.com/stblake/mathilda/blob/main/src/datetime.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
