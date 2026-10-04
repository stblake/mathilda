### Worked examples

```mathematica
In[1]:= DateList[{2026, 10, 4, 14, 30}]  (* trailing fields elide; seconds default to 0 *)
```

```mathematica
In[1]:= DateList[AbsoluteTime[{2000, 1, 1}]]  (* DateList inverts AbsoluteTime *)
```

```mathematica
In[1]:= DateList[{2022, 2, 31}]  (* out-of-range fields normalise: Feb 31 -> Mar 3 *)
```

```mathematica
In[1]:= DateList[{"9/28/2026", {"Month", "/", "Day", "/", "Year"}}]  (* explicit format elements *)
```

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
