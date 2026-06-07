# Part

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
expr[[i]] or Part[expr, i]
    gives the i-th part of expr.
expr[[-i]]
    counts from the end.
expr[[0]]
    gives the head of expr.
expr[[i, j, ...]] or Part[expr, i, j, ...]
    is equivalent to expr[[i]][[j]]..., descending into nested parts.
expr[[{i1, i2, ...}]]
    gives a list of the parts i1, i2, ... of expr (wrapped in the head of expr).
expr[[m;;n]] / expr[[m;;n;;s]]
    gives the span of parts m through n (with optional step s); ;; alone or All means all parts.

Part is treated as atomic on Integer, Real, String, Symbol, Rational[n, d], and Complex[re, im]; Part[atom, i] for i != 0 stays unevaluated.
Indices are 1-based and may be negative; out-of-range indices leave the expression unevaluated.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= {a, b, c, d}[[2]]
Out[1]= b

In[2]:= {a, b, c, d}[[-1]]
Out[2]= d

In[3]:= 123[[0]]
Out[3]= Integer
```

## Implementation notes

- Supports negative indices to count from the end (`-1` is the last element).
- `expr[[0]]` returns the `Head` of the expression. This is permitted even for atomic expressions.
- Mapping `All` across a dimension allows column extraction from matrices.

**Attributes:** `NHoldRest`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
