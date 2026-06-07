# FreeQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FreeQ[expr, form]
    yields True if no subexpression of expr matches form, False otherwise.
FreeQ[expr, form, levelspec]
    restricts the search to parts of expr at the levels specified by
    levelspec.
FreeQ[form]
    is the operator form: FreeQ[form][expr] == FreeQ[expr, form].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FreeQ[{1, 2, 4, 1, 0}, 0]
Out[1]= False

In[2]:= FreeQ[{a, b, b, a, a, a}, _Integer]
Out[2]= True

In[3]:= {f[3 x, x], f[a x, x], f[(1 + x) x, x]}
Out[3]= {3 f[x, x], a f[x, x], f[x (1 + x), x]}
```

## Implementation notes

`builtin_freeq` (`src/funcprog.c`) returns `True` iff no subexpression of the first argument matches the pattern `form`. `freeq_at_level` walks the tree depth-first, calling the pattern matcher `match` at each level permitted by the level spec (default `{0, Infinity}`), and short-circuits to `False` on the first match. `Rational`/`Complex` are treated as atomic; `Heads -> False` excludes function heads from the search.

- `Protected`.
- By default, explores levels `{0, Infinity}` and option `Heads -> True` is enabled.
- `form` can be a structural pattern.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
