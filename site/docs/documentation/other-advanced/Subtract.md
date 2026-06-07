# Subtract

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
x - y or Subtract[x, y] represents x - y; rewritten by the evaluator
to Plus[x, Times[-1, y]] so it inherits Plus's flattening and ordering.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_subtract` is a thin two-argument rewrite: `a - b` becomes `Plus[a, Times[-1, b]]`. It does no arithmetic itself — the returned `Plus`/`Times` tree is canonicalised and folded by the evaluator's `Plus`/`Times` machinery. Non-binary calls return `NULL`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/arithmetic.c`](https://github.com/stblake/mathilda/blob/main/src/arithmetic.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= 7 - 3
Out[1]= 4

In[2]:= Subtract[7, 3]
Out[2]= 4

In[3]:= a - b - c
Out[3]= a - b - c
```

### Notes

`x - y` is rewritten to `Plus[x, Times[-1, y]]`, so subtraction inherits Plus's flattening and canonical ordering rather than existing as a distinct head in the result.
