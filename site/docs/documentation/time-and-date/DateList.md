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

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= DateList[{2026, 10, 4}]
Out[1]= {2026, 10, 4, 0, 0, 0.0}

In[2]:= DateList[0]
Out[2]= {1900, 1, 1, 0, 0, 0.0}

In[3]:= DateList[{2022, 2, 31}]
Out[3]= {2022, 3, 3, 0, 0, 0.0}

In[4]:= DateList["28 Sep, 2026"]
Out[4]= {2026, 9, 28, 0, 0, 0.0}
```

### Applications (4)

Trailing fields elide; seconds default to 0

```mathematica
In[5]:= DateList[{2026, 10, 4, 14, 30}]
Out[5]= {2026, 10, 4, 14, 30, 0.0}
```

DateList inverts AbsoluteTime

```mathematica
In[6]:= DateList[AbsoluteTime[{2000, 1, 1}]]
Out[6]= {2000, 1, 1, 0, 0, 0.0}
```

Out-of-range fields normalise: Feb 31 -> Mar 3

```mathematica
In[7]:= DateList[{2022, 2, 31}]
Out[7]= {2022, 3, 3, 0, 0, 0.0}
```

Explicit format elements

```mathematica
In[8]:= DateList[{"9/28/2026", {"Month", "/", "Day", "/", "Year"}}]
Out[8]= {2026, 9, 28, 0, 0, 0.0}
```

## Implementation notes

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

- `Protected`. Returns five integers and a real second, so a fractional or current-time result carries sub-second precision (e.g. `{2026, 9, 28, 12, 24, 10.196}`).
- The inverse of `AbsoluteTime`: `DateList[AbsoluteTime[spec]]` and `DateList[spec]` agree, and out-of-range/fractional fields normalize identically (`DateList[{2022,0}] == {2021,12,1,0,0,0.}`, `DateList[{2022,1,0}] == {2021,12,31,0,0,0.}`, `DateList[{2022,3,15.5}] == {2022,3,15,12,0,0.}`).
- Year and month must be integer-valued; a non-integer month gives `DateList::arg` and leaves the input unevaluated. Day, hour, minute, and second may be noninteger.
- **Strings**: tokens are separated by any non-alphanumeric characters. `DateList["28 Sep, 2026"] == {2026,9,28,0,0,0.}`. A purely numeric string with no four-digit year is order-ambiguous: the US `Month/Day/Year` order is used and `DateList::ambig` is issued (`DateList["05/10/1"] == {2001,5,10,0,0,0.}`, two-digit years mapping to `2000+yy`); a leading four-digit token reads as ISO `Year/Month/Day`.
- **Format elements**: `"Year"`, `"YearShort"`, `"Quarter"`, `"Month"`, `"MonthName"`, `"Day"`, `"DayName"`, `"Hour"`, `"Hour12"`, `"AMPM"`, `"Minute"`, `"Second"`, `"Millisecond"`. They are read in the order given; any non-element string between them is treated as a separator. Unfilled fields default to `{current year, 1, 1, 0, 0, 0}`. E.g. `DateList[{"09/28/26",{"Day","Month","YearShort"}}] == {2028,4,9,0,0,0.}` (month `28` reduces to April 2028), `DateList[{"9/28/2026",{"Month","/","Day","/","Year"}}] == {2026,9,28,0,0,0.}`, `DateList[{"2/15",{"Month","Day"}}]` fills the current year.
- Performs no corrections for time zones, daylight saving time, or leap seconds.

**Attributes:** `Protected`.

## References

**See also:** [AbsoluteTime](../../time-and-date/AbsoluteTime/)

- Source: [`src/datetime.c`](https://github.com/stblake/mathilda/blob/main/src/datetime.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)
- Tests: [`tests/test_datetime.c`](https://github.com/stblake/mathilda/blob/main/tests/test_datetime.c)

## Notes & additional examples

### Notes

`DateList` gives a broken-down date as `{year, month, day, hour, minute,
second}` — five integers and a `Real` second, so a fractional or current-time
result keeps sub-second precision. Every form reduces to the same calendar
pipeline as `AbsoluteTime`, so `DateList[AbsoluteTime[spec]]` round-trips and
out-of-range or fractional fields normalise identically
(`DateList[{2022, 0}] == {2021, 12, 1, 0, 0, 0.}`,
`DateList[{2022, 3, 15.5}] == {2022, 3, 15, 12, 0, 0.}`).

Year and month must be integer-valued (the lengths of years and months vary); a
non-integer month gives `DateList::arg`. Free-form strings are tokenised on any
non-alphanumeric run — a purely numeric string with no four-digit year is
`Month/Day/Year` by default and warns `DateList::ambig`, while a leading
four-digit token reads as ISO `Year/Month/Day`. No time-zone, daylight-saving,
or leap-second corrections are applied.
