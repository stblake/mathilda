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
In[1]:= Accumulate[<|"a" -> 1, "b" -> 2, "c" -> 3|>]
Out[1]= <|"a" -> 1, "b" -> 3, "c" -> 6|>

In[2]:= Differences[<|"a" -> 1, "b" -> 4, "c" -> 9|>]
Out[2]= <|"b" -> 3, "c" -> 5|>

In[3]:= Ratios[<|"a" -> 1, "b" -> 2, "c" -> 6|>]
Out[3]= <|"b" -> 2, "c" -> 3|>

In[4]:= FoldList[Times, <|"a" -> 2, "b" -> 3, "c" -> 4|>]
Out[4]= <|"a" -> 2, "b" -> 6, "c" -> 24|>
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

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= FoldList[Plus, 0, {1, 2, 3, 4}]
Out[1]= {0, 1, 3, 6, 10}
```

```mathematica
In[1]:= FoldList[Times, {1, 2, 3, 4, 5}]
Out[1]= {1, 2, 6, 24, 120}
```

```mathematica
In[1]:= FoldList[Max, {3, 1, 4, 1, 5, 9, 2, 6}]
Out[1]= {3, 3, 4, 4, 5, 9, 9, 9}
```

```mathematica
In[1]:= FoldList[(#1 + #2)/2 &, 0, {1, 1, 1, 1}]
Out[1]= {0, 1/2, 3/4, 7/8, 15/16}
```

### Notes

`FoldList[f, x, list]` returns every intermediate accumulator, giving a
length-`n+1` list. `FoldList[Plus, 0, l]` produces cumulative sums and
`FoldList[Times, l]` the running partial factorials. `FoldList[Max, l]`
yields the running maximum (a prefix-scan idiom), and the exact-arithmetic
midpoint iteration `(#1 + #2)/2 &` converges on the dyadic rationals
`1 - 2^-k`, kept exact as `Rational`s rather than floats.
