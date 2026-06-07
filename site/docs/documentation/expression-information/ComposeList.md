# ComposeList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ComposeList[{f1, f2, ...}, x]
    generates a list of the form {x, f1[x], f2[f1[x]], ...}.

ComposeList applies its functions innermost-first and accumulates
the intermediate results. The output list has one more element than
the input list of functions. Function applications are evaluated
in the normal way after construction:
    ComposeList[{a, b, c}, x]  ->  {x, a[x], b[a[x]], c[b[a[x]]]}.

ComposeList has the attribute Protected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ComposeList[{a, b, c, d}, x]
Out[1]= {x, a[x], b[a[x]], c[b[a[x]]], d[c[b[a[x]]]]}

In[2]:= ComposeList[{f, g}[[{1, 2, 1, 1, 2}]], x]
Out[2]= {x, f[x], g[f[x]], f[g[f[x]]], f[f[g[f[x]]]], g[f[f[g[f[x]]]]]}

In[3]:= ComposeList[{1 - # &, 1/# &}[[{2, 2, 1, 2, 2, 1}]], x]
Out[3]= {x, 1/x, x, 1 - x, 1/(1 - x), 1 - x, x}
```

## Implementation notes

`builtin_compose_list` (`src/core.c`) takes `{f1,...,fn}` and `x` and builds the length-`n+1` list `{x, f1[x], f2[f1[x]], ...}` by constructing each symbolic application `fi[prev]`; the outer evaluator then reduces those applications to fixed point. Returns `NULL` if the first argument is not a `List`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
