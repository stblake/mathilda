# $TimeUnit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

$TimeUnit gives the minimum time interval in seconds recorded on the computer system.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= $TimeUnit
Out[1]= 1e-06

In[2]:= Head[$TimeUnit]
Out[2]= Real
```

### Applications (2)

The timer resolution, 1.*10^-6 on a standard POSIX host

```mathematica
In[3]:= $TimeUnit
Out[3]= 1e-06
```

Always a small positive machine real

```mathematica
In[4]:= $TimeUnit > 0
Out[4]= True
```

## Options & behaviour

The value is machine-dependent, but on a standard POSIX host `CLOCKS_PER_SEC` is
`1000000`, so `$TimeUnit` is `1.*10^-6`.

## Implementation notes

**Algorithm.** `$TimeUnit` is not a function but a read-only system constant,
registered in `system_constants_init` (`src/core.c`) by
`register_system_constant("$TimeUnit", expr_new_real(1.0 / (double)CLOCKS_PER_SEC))`.
Its value is the reciprocal of the C library's `CLOCKS_PER_SEC`, i.e. the
resolution of the `clock()`-based CPU timers that `Timing` and `TimeUsed` read,
and the granularity `Pause` documents itself against. On a standard POSIX host
`CLOCKS_PER_SEC` is `1000000`, so `$TimeUnit` is `1.*10^-6`.

**Data structures.** A single `EXPR_REAL` bound as the symbol's `OwnValue`;
nothing is recomputed per reference. `register_system_constant` also sets
`ATTR_PROTECTED`, so the binding cannot be reassigned.

**Complexity / limits.** O(1) lookup. The value is a machine `double` fixed at
start-up from a compile-time constant; it is a scalar system parameter, not a
numeric kernel, so there is no packed/`Compile[]` surface to provide.

- `Protected` (read-only system constant).
- A real equal to the resolution of the `clock()`-based timers
  (`1 / CLOCKS_PER_SEC`), the granularity `Pause` documents itself against.

**Attributes:** `Protected`.

## References

**See also:** [Pause](../../time-and-date/Pause/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/time-and-date.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/time-and-date.md)

## Notes & additional examples

### Notes

`$TimeUnit` is a read-only system constant: the minimum time interval the
computer records, equal to `1/CLOCKS_PER_SEC`. It is the resolution of the
`clock()`-based CPU timers behind `Timing` and `TimeUsed`, and the granularity
`Pause` is accurate down to. Its value is machine-dependent, but on standard
systems `CLOCKS_PER_SEC` is `1000000`, so `$TimeUnit` is `1.*10^-6`. Being
`Protected`, it cannot be reassigned.
