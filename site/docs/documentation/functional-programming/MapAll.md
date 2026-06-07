# MapAll

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
f //@ expr or MapAll[f, expr]
    applies f to every subexpression in expr (equivalent to
    Map[f, expr, {0, Infinity}]).  Atomic leaves are wrapped too.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_map_all` is a thin wrapper around the same `map_at_level` traversal used
by `Map`, but with the fixed level-spec `{0, Infinity}` (`min=0`,
`max=1000000`, `heads=false`), i.e. `MapAll[f, expr]` ≡ `Map[f, expr, {0,
Infinity}]`. The bottom-up recursion rebuilds every `EXPR_FUNCTION` from its
mapped children and then wraps each node (including the whole expression at level
0) in `f[...]`, calling `evaluate()` so `f`'s attributes apply. A trailing
`Heads -> True` option is honoured via `parse_options`.

**Attributes:** `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= MapAll[f, {a, {b, c}}]
Out[1]= f[{f[a], f[{f[b], f[c]}]}]

In[2]:= f //@ {a, b}
Out[2]= f[{f[a], f[b]}]

In[3]:= MapAll[g, 1 + x]
Out[3]= g[g[1] + g[x]]

In[4]:= Map[f, {a, {b, c}}]
Out[4]= {f[a], f[{b, c}]}
```

### Notes

`MapAll` applies `f` to every subexpression including atomic leaves, equivalent to `Map[f, expr, {0, Infinity}]`; its operator form is `f //@ expr`. Unlike `Map`, which only touches the first level (compare In[1] vs In[4]), `MapAll` reaches all levels and wraps the whole expression too.
