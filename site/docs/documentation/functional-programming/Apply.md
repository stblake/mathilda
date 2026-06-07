# Apply

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
f @@ expr or Apply[f, expr]
    replaces the head of expr with f.
Apply[f, expr, levelspec]
    performs the head replacement at the parts of expr specified by
    levelspec; the default levelspec is {0} (top level only).
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Apply[Plus, {1, 2, 3, 4}]
Out[1]= 10
```

```mathematica
In[1]:= f @@ {a, b, c}
Out[1]= f[a, b, c]
```

```mathematica
In[1]:= Apply[List, a + b + c]
Out[1]= {a, b, c}
```

```mathematica
In[1]:= Apply[f, {{a, b}, {c, d}}, {1}]
Out[1]= {f[a, b], f[c, d]}
```

### Notes

`f @@ expr` is the shorthand for `Apply[f, expr]`: it replaces the head of `expr`
with `f`. `Apply[Plus, list]` is the standard idiom for summing a list, since
the list's `List` head is swapped for `Plus`. Because addition is stored as
`Plus[...]`, `Apply[List, a + b + c]` recovers the summands. A level
specification like `{1}` applies the head replacement to each element at that
level instead of the whole expression.
