# Order

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Order[e1, e2] gives 1 if e1 is before e2 in canonical order, -1 if e1 is after e2, and 0 if e1 is identical to e2.`**

<details>
<summary>Notes</summary>

Order compares structurally (the same canonical order as Sort), not by numerical value, and is compilable.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= {Order[a, a], Order[a, b], Order[b, a]}
Out[1]= {0, 1, -1}

In[2]:= {Order[6, Pi], Order[6, N[Pi]]}
Out[2]= {1, -1}

In[3]:= Order @@@ Tuples[{0, 1, 2}, 2]
Out[3]= {0, 1, 1, -1, 0, 1, -1, -1, 0}
```

### Applications (5)

A sorts before b: +1

```mathematica
In[4]:= Order[a, b]
Out[4]= 1
```

Reversed: -1

```mathematica
In[5]:= Order[b, a]
Out[5]= -1
```

Identical: 0

```mathematica
In[6]:= Order[a, a]
Out[6]= 0
```

Structural: the Integer 6 sorts before the symbol Pi

```mathematica
In[7]:= Order[6, Pi]
Out[7]= 1
```

Two numeric atoms, now compared by value

```mathematica
In[8]:= Order[6, N[Pi]]
Out[8]= -1
```

## Implementation notes

**Algorithm.** `builtin_order` is the user-facing surface of `expr_compare`, the
canonical comparator every sorting routine (`Sort`, `SortBy`, `OrderedQ`,
`Ordering`, ...) is built on. `Order[e1, e2]` requires exactly two arguments and
returns `1` if `e1` is before `e2` in canonical order, `-1` if after, `0` if
identical. `expr_compare` returns negative when `e1` sorts first, so the sign is
simply inverted. Comparison is **structural**, not by numerical value:
`Order[6, Pi]` is `1` (the `Integer` `6` sorts before the symbol `Pi`) whereas
`Order[6, N[Pi]]` is `-1` (two numeric atoms, compared by value) — the canonical
order ranks reals by value, then strings, then symbols, then expressions by
length/head/parts.

**Data structures.** None of its own — one `expr_compare` call over the two
borrowed argument trees, returning an `EXPR_INTEGER` in `{1, 0, -1}`.

**Complexity / limits.** `O(min size)` of the two trees in the worst case (the
comparator short-circuits on the first difference). It **is** compilable:
`CompileDiagnostics[{{x, _Real}}, Order[x, 1]]` reports `Compiled -> True` with
`ResultType -> Integer`, lowering over machine numbers to `Sign[e2 - e1]`
(matching the interpreter's integer head); complex or array arguments fall back
to the interpreter. No packed/NDArray path — it is a scalar two-argument
comparator, not an element-wise map. `Protected`.

- `Protected`.
- Uses the same internal canonical comparison (`expr_compare`) as `Sort` and `OrderedQ` — see the canonical-order rules under `Sort` below.
- Compares **structurally**, not by numerical value: `Order[6, Pi]` is `1` (the Integer `6` sorts before the symbol `Pi`), whereas `Order[6, N[Pi]]` is `-1` (two numeric atoms, compared by value).
- Compilable inside `Compile[]` and auto-compiled by `Plot`/`Table`/`NIntegrate`: over machine numbers it lowers to `Sign[e2 - e1]`, returning the Integer `{1, 0, -1}` (matching the interpreter's head). Complex/array arguments fall back to the interpreter.
- Requires exactly two arguments; otherwise it stays unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Sort](../../data-structures/Sort/), [OrderedQ](../../structural-manipulation/OrderedQ/), [Pi](../../mathematical-constants/Pi/), [Plot](../../graphics/Plot/), [Table](../../lists-and-iteration/Table/), [NIntegrate](../../numerical-calculus/NIntegrate/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_assoc_read.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_read.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_sort.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sort.c)

## Notes & additional examples

### Notes

`Order[e1, e2]` is the user-facing surface of the canonical comparator every
sorting routine (`Sort`, `SortBy`, `OrderedQ`, `Ordering`) is built on: `1` if
`e1` is before `e2`, `-1` if after, `0` if identical. The comparison is
**structural**, not by numerical value — `Order[6, Pi]` is `1` because an
`Integer` sorts before a symbol, whereas `Order[6, N[Pi]]` is `-1` because two
numeric atoms compare by value. It needs exactly two arguments (otherwise it
stays unevaluated) and is compilable: over machine numbers it lowers to
`Sign[e2 - e1]`, returning the same `{1, 0, -1}` integer.
