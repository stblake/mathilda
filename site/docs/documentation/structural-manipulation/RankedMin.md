# RankedMin

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RankedMin[list, n]`**

gives the n-th smallest element of list.

**`RankedMin[list, -n]`**

gives the n-th largest element of list. RankedMin\[list, 1\] is Min\[list\] and RankedMin\[list, -1\] is Max\[list\]. Yields a definite result when every element is a real number; +-Infinity are ordered as +-infinity. Has a packed-array fast path and is compilable.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= RankedMin[{12, 13, 11}, 2]
Out[1]= 12

In[2]:= RankedMin[{Pi, Sqrt[2], E, 3}, 3]
Out[2]= 3

In[3]:= RankedMax[{2.5, E, 12, 15, 485}, -2]
Out[3]= E

In[4]:= RankedMax[{Infinity, 5, Infinity, -Infinity}, 2]
Out[4]= Infinity
```

### Applications (3)

The 2nd smallest element

```mathematica
In[5]:= RankedMin[{12, 13, 11}, 2]
Out[5]= 12
```

Symbolic reals order by value

```mathematica
In[6]:= RankedMin[{Pi, Sqrt[2], E, 3}, 3]
Out[6]= 3
```

The 1st largest is Max

```mathematica
In[7]:= RankedMin[{12, 13, 11}, -1]
Out[7]= 13
```

## Implementation notes

**Algorithm.** `builtin_ranked_min` selects the `n`-th smallest element of a list
— `RankedMin[list, n]`, with `RankedMin[list, -n]` the `n`-th largest, so
`RankedMin[list, 1]` is `Min[list]` and `RankedMin[list, -1]` is `Max[list]`. It
requires an integer, nonzero `n`, converts it to a 1-based ascending rank `r`
(`r = n` if positive, `m + n + 1` if negative; an out-of-range `r` or empty list
leaves the call unevaluated), and runs `ranked_select`. One scan chooses the
comparison path: if every element is a real numeric atom it compares exactly with
`expr_compare` (BigInt/Rational-safe), otherwise it builds a `double` key per
element with `ranked_numeric_key` (`±HUGE_VAL` for `±Infinity`, the value of a
real atom, else the machine-precision numericalisation of a symbolic real such as
`Pi`, `E`, `Sqrt[2]`) — and a non-real element (free symbol, non-real complex)
makes it decline, since a definite result needs every element to be a real
number. `ranked_select_idx` then quickselects the `r`-th position with a Hoare
partition and middle pivot (`O(m)` average), the original index breaking value
ties stably, and the exact element at that position is returned.

**Data structures.** A `size_t* idx` index array quickselected in place; an
optional `double* keys` for the approximate path; a `RankedCtx {Expr** elem;
const double* keys;}` read by `ranked_cmp`. Only the `r`-th slot is guaranteed
settled on return. The result is an `expr_copy` of the selected element in its
exact form.

**Complexity / limits.** `O(m)` average quickselect. An `NDArray` argument takes
the buffer order-statistic path `ndred_ranked_min` (`int64` exactly, reals via an
`O(m)` quickselect) — `RankedMin` is on `pack.c`'s `AWARE` and `INT64_OK` lists.
It is also compilable:
`CompileDiagnostics[{{v, _Real, 1}}, RankedMin[v, 2]]` reports `Compiled -> True`.
`Protected`.

- `Protected`.
- `RankedMax[list, k]` is `RankedMin[list, -k]`.
- `RankedMin[list, 1]` is `Min[list]`; `RankedMin[list, -1]` is `Max[list]`.
- Yields a definite result whenever every element is a real number, including
  symbolic real constants (`Pi`, `E`, `Sqrt[2]`, `Pi + E`), which order by value;
  `Infinity`/`-Infinity` rank as `±∞`. Returns the element in its exact form.
- Exact for arbitrary-precision integers and rationals; a symbolic non-real
  element (a free symbol or a non-real complex), an empty list, or `|n|` out of
  range leaves the call unevaluated.
- Packed-array fast path (int64 exact, real via O(*n*) quickselect) and a
  `Compile[]` lowering, so `RankedMin[v, k]`/`RankedMax[v, k]` compile and
  auto-compile.

**Attributes:** `Protected`.

## References

**See also:** [RankedMax](../../structural-manipulation/RankedMax/), [Min](../../data-structures/Min/), [Max](../../data-structures/Max/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_ranked.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ranked.c)

## Notes & additional examples

### Notes

`RankedMin[list, n]` is the `n`-th smallest element — the order statistic between
`Min` and `Max`. A negative index counts from the top, so `RankedMin[list, -n]`
is the `n`-th largest, `RankedMin[list, 1]` is `Min[list]`, and
`RankedMin[list, -1]` is `Max[list]`. It returns a definite result whenever every
element is a real number, including symbolic real constants (`Pi`, `E`,
`Sqrt[2]`), which order by value, with `±Infinity` ranking as `±∞`; the element
comes back in its exact form. A non-real element, an empty list, or `|n|` out of
range leaves the call unevaluated. Selection is an `O(n)` quickselect (int64
exact), with a packed-array fast path and a `Compile[]` lowering.
