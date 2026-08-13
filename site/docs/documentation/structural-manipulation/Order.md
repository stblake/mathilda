# Order

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Order[e1, e2] gives 1 if e1 is before e2 in canonical order, -1 if e1 is after e2, and 0 if e1 is identical to e2.`**

<details>
<summary>Notes</summary>

Order compares structurally (the same canonical order as Sort), not by numerical value, and is compilable.

</details>

## Examples (3)

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

## Implementation notes

- `Protected`.
- Uses the same internal canonical comparison (`expr_compare`) as `Sort` and `OrderedQ` — see the canonical-order rules under `Sort` below.
- Compares **structurally**, not by numerical value: `Order[6, Pi]` is `1` (the Integer `6` sorts before the symbol `Pi`), whereas `Order[6, N[Pi]]` is `-1` (two numeric atoms, compared by value).
- Compilable inside `Compile[]` and auto-compiled by `Plot`/`Table`/`NIntegrate`: over machine numbers it lowers to `Sign[e2 - e1]`, returning the Integer `{1, 0, -1}` (matching the interpreter's head). Complex/array arguments fall back to the interpreter.
- Requires exactly two arguments; otherwise it stays unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [Sort](../../data-structures/Sort/), [OrderedQ](../../structural-manipulation/OrderedQ/), [Pi](../../mathematical-constants/Pi/), [Plot](../../graphics/Plot/), [Table](../../lists-and-iteration/Table/), [NIntegrate](../../numerical-calculus/NIntegrate/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_sort.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sort.c)
