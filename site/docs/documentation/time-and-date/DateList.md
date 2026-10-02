# DateList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DateList[]`**

gives the current local date and time in the form {y, m, d, h, m, s}.

**`DateList[date]`**

gives a date list corresponding to the given date specification.

<details>
<summary>Notes</summary>

The supported date specifications are: {y, m, d, h, m, s}    DateList specification time            AbsoluteTime specification (a number of seconds since 1900) "string"        DateString specification {"string", {e1, ...}}    date string parsed with the given format elements In {y, m, ...} entries may be elided from the right: {y} is {y,1,1,0,0,0}, {y,m} is {y,m,1,0,0,0}, and so on. Values of m, d, h, m, s outside their normal ranges are reduced (m=0 is the previous December, d=0 the last day of the previous month); d, h, m, s may be noninteger, but the year and month must be integers. The format elements are "Year", "YearShort", "Quarter", "Month", "MonthName", "Day", "DayName", "Hour", "Hour12", "AMPM", "Minute", "Second", and "Millisecond". Elements are read from the string in the order given, separated by any non-alphanumeric characters (or by explicit separator strings placed between them). Fields not filled default to {current year, 1, 1, 0, 0, 0}. DateList uses the local date and time with no correction for time zones, daylight saving time, or leap seconds.

</details>

## Examples

_No verified examples yet for this function._

## Implementation notes

- `Protected`. Returns five integers and a real second, so a fractional or current-time result carries sub-second precision (e.g. `{2026, 9, 28, 12, 24, 10.196}`).
- The inverse of `AbsoluteTime`: `DateList[AbsoluteTime[spec]]` and `DateList[spec]` agree, and out-of-range/fractional fields normalize identically (`DateList[{2022,0}] == {2021,12,1,0,0,0.}`, `DateList[{2022,1,0}] == {2021,12,31,0,0,0.}`, `DateList[{2022,3,15.5}] == {2022,3,15,12,0,0.}`).
- Year and month must be integer-valued; a non-integer month gives `DateList::arg` and leaves the input unevaluated. Day, hour, minute, and second may be noninteger.
- **Strings**: tokens are separated by any non-alphanumeric characters. `DateList["28 Sep, 2026"] == {2026,9,28,0,0,0.}`. A purely numeric string with no four-digit year is order-ambiguous: the US `Month/Day/Year` order is used and `DateList::ambig` is issued (`DateList["05/10/1"] == {2001,5,10,0,0,0.}`, two-digit years mapping to `2000+yy`); a leading four-digit token reads as ISO `Year/Month/Day`.
- **Format elements**: `"Year"`, `"YearShort"`, `"Quarter"`, `"Month"`, `"MonthName"`, `"Day"`, `"DayName"`, `"Hour"`, `"Hour12"`, `"AMPM"`, `"Minute"`, `"Second"`, `"Millisecond"`. They are read in the order given; any non-element string between them is treated as a separator. Unfilled fields default to `{current year, 1, 1, 0, 0, 0}`. E.g. `DateList[{"09/28/26",{"Day","Month","YearShort"}}] == {2028,4,9,0,0,0.}` (month `28` reduces to April 2028), `DateList[{"9/28/2026",{"Month","/","Day","/","Year"}}] == {2026,9,28,0,0,0.}`, `DateList[{"2/15",{"Month","Day"}}]` fills the current year.
- Performs no corrections for time zones, daylight saving time, or leap seconds.

**Attributes:** `Protected`.

## References

**See also:** [AbsoluteTime](../../time-and-date/AbsoluteTime/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_datetime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_datetime.c)
