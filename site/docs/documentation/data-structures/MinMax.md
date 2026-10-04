# MinMax

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MinMax[list]`**

Gives {Min\[list\], Max\[list\]}. Over an association, uses the values.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= MinMax[<|"a" -> 3, "b" -> 1, "c" -> 9|>]
Out[1]= {1, 9}
```

### Applications (2)

{smallest, largest} in one shot

```mathematica
In[2]:= MinMax[{3, 1, 4, 1, 5, 9, 2, 6}]
Out[2]= {1, 9}
```

Over an association, uses the values

```mathematica
In[3]:= MinMax[<|a -> 3, b -> 1, c -> 4|>]
Out[3]= {1, 4}
```

## Implementation notes

**Algorithm.** `builtin_minmax` returns `{Min[arg], Max[arg]}`, built by evaluating a
`Min` call and a `Max` call on (refcount-shared copies of) the argument. Delegating
keeps every numeric subtlety — bignums, reals, symbolic extrema, the empty-list
`Infinity`/`-Infinity` — in exactly one place rather than reimplemented here. A
packed buffer, an `NDArray` and an association are all accepted; over an association
the values are used (as `Min`/`Max` already do).

**Data structures.** No bespoke storage: `expr_copy` is a refcount bump, so the same
buffer is handed to both `Min` and `Max`, and each takes its own buffer reduction
(`ndred_min` / `ndred_max`).

**Complexity / limits.** O(n) — effectively two single-pass reductions over the same
data. `MinMax` is on `pack.c`'s `AWARE` list, so a packed/`NDArray` argument stays on
the buffer; without that both halves would materialise 10⁶ boxed nodes (a measured
~735× regression against the buffer path). For a non-list/array/association argument
it declines and the call is left unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Min](../../data-structures/Min/), [Max](../../data-structures/Max/)

- Source: [`src/list/minmax.c`](https://github.com/stblake/mathilda/blob/main/src/list/minmax.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)

## Notes & additional examples

### Notes

`MinMax[list]` is `{Min[list], Max[list]}` computed without writing the two calls by
hand. It threads onto the packed-array fast path, so a large numeric vector stays on
the buffer and both extrema come from one data pass each rather than from 10⁶ boxed
elements.
