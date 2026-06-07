# FoldList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FoldList[f, x, list]
    gives {x, f[x, list[[1]]], f[f[x, list[[1]]], list[[2]]], ...}.
FoldList[f, list]
    gives {list[[1]], f[list[[1]], list[[2]]], ...}.

For a length-n list, FoldList generates a list of length n+1. The
head of list is preserved in the output:
    FoldList[f, x, p[a, b]] -> p[x, f[x, a], f[f[x, a], b]].
FoldList[f, {}] returns an empty list {}. f may be a symbol or a
pure function; each intermediate application is evaluated before the
next one.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FoldList[f, x, {a, b, c, d}]
Out[1]= {x, f[x, a], f[f[x, a], b], f[f[f[x, a], b], c], f[f[f[f[x, a], b], c], d]}

In[2]:= FoldList[f, {a, b, c, d}]
Out[2]= {a, f[a, b], f[f[a, b], c], f[f[f[a, b], c], d]}

In[3]:= FoldList[Plus, 0, Range[5]]
Out[3]= {0, 1, 3, 6, 10, 15}

In[4]:= FoldList[f, x, p[a, b, c, d]]
Out[4]= p[x, f[x, a], f[f[x, a], b], f[f[f[x, a], b], c], f[f[f[f[x, a], b], c], d]]

In[5]:= FoldList[g[#2, #1] &, x, {a, b, c, d}]
Out[5]= {x, g[a, x], g[b, g[a, x]], g[c, g[b, g[a, x]]], g[d, g[c, g[b, g[a, x]]]]}

In[6]:= FoldList[Times, 1, Range[10]]
Out[6]= {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800}

In[7]:= FoldList[1/(#2 + #1) &, x, Reverse[{a, b, c}]]
Out[7]= {x, 1/(c + x), 1/(b + 1/(c + x)), 1/(a + 1/(b + 1/(c + x)))}

In[8]:= FoldList[10 #1 + #2 &, 0, {4, 5, 1, 6, 7, 8}]
Out[8]= {0, 4, 45, 451, 4516, 45167, 451678}
```

## Implementation notes

`builtin_foldlist` is `fold_impl(res, true)`: it returns every intermediate
accumulator, `FoldList[f, x, {a,b,c}]` → `{x, f[x,a], f[f[x,a],b],
f[f[f[x,a],b],c]}`. The two-argument form seeds with the list's first element;
`FoldList[f, {}]` returns an empty list under the input list's head. The seed is
pushed into an `ExprBuf`, and `iter_run` with `fold_step` consumes the remaining
elements left-to-right via `apply_binary(f, acc, elem)` (`eval_and_free`). With
`as_list=true`, `ebuf_finalize` wraps the full history under the **input list's
head** (preserved via `expr_copy(list_head)`). Shares all machinery with `Fold`.

- `Protected`.
- With a length-`n` list, `FoldList` returns a list of length `n + 1` (or `n` in the no-seed form).
- The head of the third argument is preserved in the output: `FoldList[f, x, p[a, b]]` gives `p[x, f[x, a], f[f[x, a], b]]`.
- `FoldList[f, {}]` returns `{}` (an empty list with the input head); `FoldList[f, x, {}]` returns `{x}`.
- `Fold[f, x, list]` is equivalent to `Last[FoldList[f, x, list]]`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
