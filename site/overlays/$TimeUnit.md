### Worked examples

```mathematica
In[1]:= $TimeUnit  (* the timer resolution, 1.*10^-6 on a standard POSIX host *)
```

```mathematica
In[1]:= $TimeUnit > 0  (* always a small positive machine real *)
```

### Notes

`$TimeUnit` is a read-only system constant: the minimum time interval the
computer records, equal to `1/CLOCKS_PER_SEC`. It is the resolution of the
`clock()`-based CPU timers behind `Timing` and `TimeUsed`, and the granularity
`Pause` is accurate down to. Its value is machine-dependent, but on standard
systems `CLOCKS_PER_SEC` is `1000000`, so `$TimeUnit` is `1.*10^-6`. Being
`Protected`, it cannot be reassigned.
