# MapAll

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

f //@ expr or MapAll\[f, expr\] applies f to every subexpression in expr (equivalent to Map\[f, expr, {0, Infinity}\]).  Atomic leaves are wrapped too.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (6)

```mathematica
In[1]:= MapAll[f, {a, {b, c}}]
Out[1]= f[{f[a], f[{f[b], f[c]}]}]

In[2]:= f //@ {a, b}
Out[2]= f[{f[a], f[b]}]

In[3]:= MapAll[g, 1 + x]
Out[3]= g[g[1] + g[x]]

In[4]:= Map[f, {a, {b, c}}]
Out[4]= {f[a], f[{b, c}]}

In[5]:= f //@ (x^2 + y)
Out[5]= f[f[f[x]^f[2]] + f[y]]

In[6]:= MapAll[g, 1 + x^2]
Out[6]= g[g[1] + g[g[x]^g[2]]]
```

## Implementation notes

`builtin_map_all` is a thin wrapper around the same `map_at_level` traversal used
by `Map`, but with the fixed level-spec `{0, Infinity}` (`min=0`,
`max=1000000`, `heads=false`), i.e. `MapAll[f, expr]` ≡ `Map[f, expr, {0,
Infinity}]`. The bottom-up recursion rebuilds every `EXPR_FUNCTION` from its
mapped children and then wraps each node (including the whole expression at level
0) in `f[...]`, calling `evaluate()` so `f`'s attributes apply. A trailing
`Heads -> True` option is honoured via `parse_options`.

**Attributes:** `Protected`.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_map_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_map_ndarray.c)
- Tests: [`tests/test_sow_reap.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sow_reap.c)

## Notes & additional examples

### Notes

`MapAll` applies `f` to every subexpression including atomic leaves, equivalent to `Map[f, expr, {0, Infinity}]`; its operator form is `f //@ expr`. Unlike `Map`, which only touches the first level (compare In[1] vs In[4]), `MapAll` reaches all levels and wraps the whole expression too.
