# UnixTime

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`UnixTime[]`**

gives the total number of seconds since the beginning of January 1, 1970, GMT.

**`UnixTime[date]`**

gives the Unix time corresponding to the given date specification.

**`UnixTime[] gives the number of seconds elapsed since {1970, 1, 1, 0, 0, 0} GMT, not`**

<details>
<summary>Notes</summary>

The supported date specifications are: {y, m, d, h, m, s}    DateList specification time            AbsoluteTime specification (a number of seconds since 1900) "string"        DateString specification {"string", {e1, ...}}    date string parsed with the given format elements counting leap seconds, and always returns the nearest whole second as an integer. In {y, m, ...} entries may be elided from the right ({y} is {y,1,1,0,0,0}, {y,m} is {y,m,1,0,0,0}, and so on); values outside their normal ranges are reduced and the result is rounded to the nearest second. The year and month must be integers. UnixTime interprets a date the same way DateList does and applies no correction for time zones, daylight saving time, or leap seconds; it is AbsoluteTime shifted by the fixed 1900-to-1970 epoch offset of 2208988800 seconds.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= UnixTime[{1970, 1, 1, 0, 0, 0}]
Out[1]= 0

In[2]:= UnixTime[{2022, 1, 1, 0, 0, 0}]
Out[2]= 1640995200

In[3]:= UnixTime[{2022, 2, 31}]
Out[3]= 1646265600

In[4]:= AbsoluteTime[{2000, 1, 1}] - UnixTime[{2000, 1, 1}]
Out[4]= 2208988800
```

### Applications (4)

Seconds since 1970-01-01 GMT

```mathematica
In[5]:= UnixTime[{2022, 1, 1, 0, 0, 0}]
Out[5]= 1640995200
```

The fixed 1900->1970 offset

```mathematica
In[6]:= AbsoluteTime[{2000, 1, 1}] - UnixTime[{2000, 1, 1}]
Out[6]= 2208988800
```

A date string, parsed as DateList would

```mathematica
In[7]:= UnixTime["1 Jan 2000"]
Out[7]= 946684800
```

The number 0 is an AbsoluteTime (1900), so this is negative

```mathematica
In[8]:= UnixTime[0]
Out[8]= -2208988800
```

## Implementation notes

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

- `Protected`.
- `AbsoluteTime` shifted by the fixed 1900→1970 epoch offset `days_since_1900[1970,1,1] * 86400 == 2208988800`: `UnixTime[spec] == AbsoluteTime[spec] - 2208988800`, so `UnixTime[{1970,1,1,0,0,0}] == 0` and `UnixTime[{2022,1,1,0,0,0}] == 1640995200`.
- **Always returns an integer**, the nearest whole second — a fractional-second spec is rounded rather than kept as a real (a deliberate difference from `AbsoluteTime`).
- `UnixTime[]` reports the true POSIX second (GMT), matching the operating system clock (`date +%s`); the date-list and string forms are timezone-free calendar arithmetic and interpret a spec exactly as `DateList` does.
- Year and month must be integer-valued; a non-integer month gives `UnixTime::arg` and leaves the input unevaluated. Day, hour, minute, and second may be noninteger.
- Out-of-range date components are reduced to standard normalized form before rounding, e.g. `UnixTime[{2022, 2, 31}] == UnixTime[{2022, 3, 3}] == 1646265600`.
- Performs no corrections for time zones, daylight saving time, or leap seconds.

**Attributes:** `Protected`.

## References

**See also:** [DateList](../../time-and-date/DateList/), [AbsoluteTime](../../time-and-date/AbsoluteTime/)

- Source: [`src/datetime.c`](https://github.com/stblake/mathilda/blob/main/src/datetime.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_datetime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_datetime.c)

## Notes & additional examples

### Notes

`UnixTime` is `AbsoluteTime` shifted by the fixed epoch offset
`2208988800` (= `25567 * 86400`), so `UnixTime[spec] == AbsoluteTime[spec] -
2208988800` and `UnixTime[{1970, 1, 1, 0, 0, 0}] == 0`. It **always returns an
integer** — the nearest whole second — so a fractional-second spec is rounded
rather than kept as a real (a deliberate difference from `AbsoluteTime`).

`UnixTime[]` reports the true POSIX second in GMT, matching the operating-system
clock (`date +%s`). A bare number argument is read as an `AbsoluteTime` (seconds
since 1900), which is why `UnixTime[0]` is negative. The date-list and string
forms share `DateList`'s parser, so the two heads always agree on how a spec is
interpreted; no time-zone, DST, or leap-second corrections are applied.
