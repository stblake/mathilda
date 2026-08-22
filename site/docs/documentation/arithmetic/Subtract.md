# Subtract

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

x - y or Subtract\[x, y\] represents x - y; rewritten by the evaluator to Plus\[x, Times\[-1, y\]\] so it inherits Plus's flattening and ordering.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= 7 - 3
Out[1]= 4

In[2]:= Subtract[7, 3]
Out[2]= 4

In[3]:= a - b - c
Out[3]= a - b - c
```

## Implementation notes

`builtin_subtract` is a thin two-argument rewrite: `a - b` becomes `Plus[a, Times[-1, b]]`. It does no arithmetic itself — the returned `Plus`/`Times` tree is canonicalised and folded by the evaluator's `Plus`/`Times` machinery. Non-binary calls return `NULL`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/arithmetic.c`](https://github.com/stblake/mathilda/blob/main/src/arithmetic.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`x - y` is rewritten to `Plus[x, Times[-1, y]]`, so subtraction inherits Plus's flattening and canonical ordering rather than existing as a distinct head in the result.
