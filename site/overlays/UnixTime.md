### Worked examples

```mathematica
In[1]:= UnixTime[{2022, 1, 1, 0, 0, 0}]  (* seconds since 1970-01-01 GMT *)
```

```mathematica
In[1]:= AbsoluteTime[{2000, 1, 1}] - UnixTime[{2000, 1, 1}]  (* the fixed 1900->1970 offset *)
```

```mathematica
In[1]:= UnixTime["1 Jan 2000"]  (* a date string, parsed as DateList would *)
```

```mathematica
In[1]:= UnixTime[0]  (* the number 0 is an AbsoluteTime (1900), so this is negative *)
```

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
