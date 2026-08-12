# RankedMax

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RankedMax[list, n]`**

gives the n-th largest element of list.

**`RankedMax[list, -n]`**

gives the n-th smallest element of list. RankedMax\[list, n\] is RankedMin\[list, -n\]; RankedMax\[list, 1\] is Max\[list\] and RankedMax\[list, -1\] is Min\[list\]. Yields a definite result when every element is a real number. Has a packed-array fast path and is compilable.

## Examples (4)

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

## Implementation notes

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

## See also

[RankedMin](../../structural-manipulation/RankedMin/), [Min](../../data-structures/Min/), [Max](../../data-structures/Max/), [Pi](../../mathematical-constants/Pi/), [E](../../mathematical-constants/E/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_ranked.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ranked.c)
