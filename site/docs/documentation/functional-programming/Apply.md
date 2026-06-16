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

**Algorithm.** `builtin_apply` replaces the head of `expr` with `f` at the
levels given by an optional level-spec (default level `0`, i.e. just the top
head). The work is done by the recursive `apply_at_level`, which descends each
`EXPR_FUNCTION`. A node is "at" the active range when its current depth lies in
`[spec.min, spec.max]`, with negative levels measured against `get_depth(expr)`
(so `-1` etc. count from the leaves). When a node is in range its arguments are
first transformed recursively, a fresh `f[args...]` is built, and `evaluate()`
is called on it so `f`'s attributes take effect; otherwise the original head is
kept (or transformed too when `Heads -> True`).

**Level / option parsing.** The third argument is interpreted by
`parse_level_spec` (handles an integer `n`, `{n}`, `{m,n}`, `Infinity`, and
treats `Automatic`/missing as level `{0,0}`); any trailing `Heads -> True`
option is read by `parse_options`. A `Rule`-headed third argument is recognised
as an option, not a level-spec.

**Data structures.** Operates directly on the `Expr` tagged union; rebuilds
`EXPR_FUNCTION` nodes with `expr_new_function`. Leaves (non-function atoms) are
returned via `expr_copy` unchanged.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
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

```mathematica
In[1]:= Apply[Times, Range[10]]
Out[1]= 3628800
```

```mathematica
In[1]:= Apply[GCD, {84, 126, 210}]
Out[1]= 42

In[2]:= Apply[Plus, Table[1/k^2, {k, 1, 6}]]
Out[2]= 5369/3600
```

### Notes

`f @@ expr` is the shorthand for `Apply[f, expr]`: it replaces the head of `expr`
with `f`. `Apply[Plus, list]` is the standard idiom for summing a list, since
the list's `List` head is swapped for `Plus`. Because addition is stored as
`Plus[...]`, `Apply[List, a + b + c]` recovers the summands. A level
specification like `{1}` applies the head replacement to each element at that
level instead of the whole expression. Folding a variadic head over a generated
list is a common pattern: `Apply[Times, Range[10]]` is `10!`, `Apply[GCD, ...]`
collapses a list to a single greatest common divisor, and `Apply[Plus, ...]`
over exact rationals returns the exact partial sum `5369/3600`.
