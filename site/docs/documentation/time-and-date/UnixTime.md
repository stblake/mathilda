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

## Examples

_No verified examples yet for this function._

## Implementation notes

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

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_datetime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_datetime.c)
