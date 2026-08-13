# Slot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

# or Slot\[n\] represents the n-th argument of a pure function.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= (#1 + #2 &)[3, 4]
Out[1]= 7

In[2]:= Map[#^2 &, {1, 2, 3, 4}]
Out[2]= {1, 4, 9, 16}

In[3]:= FullForm[#2]
Out[3]= Slot[2]
```

## Implementation notes

`Slot[n]` (`#`, `#n`) is an inert marker — `builtin_slot` always returns `NULL` so the node stays unevaluated on its own. Substitution happens only when a pure `Function` is applied: `substitute_slots` walks the function body and replaces each `Slot[n]` with a copy of the `n`-th argument (when `1 ≤ n ≤ arg_count`), stopping recursion at any nested `Function` so inner pure functions are not captured by the outer slots. `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## References

- Source: [`src/purefunc.c`](https://github.com/stblake/mathilda/blob/main/src/purefunc.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_minimalpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_minimalpolynomial.c)
- Tests: [`tests/test_numloop.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numloop.c)
- Tests: [`tests/test_patterns.c`](https://github.com/stblake/mathilda/blob/main/tests/test_patterns.c)

## Notes & additional examples

### Notes

`#` (or `#1`) is `Slot[1]`, the first argument of the enclosing pure function (`&`); `#n` refers to the n-th argument. A bare `#` inside `Map` (Out[2]) receives each list element in turn.
