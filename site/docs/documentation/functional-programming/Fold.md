# Fold

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Fold[f, x, list]
    gives the last element of FoldList[f, x, list]:
    f[...f[f[f[x, list[[1]]], list[[2]]], list[[3]]]..., list[[n]]].
Fold[f, list]
    is equivalent to Fold[f, First[list], Rest[list]].

The head of list need not be List. Fold[f, x, {}] returns x, and
Fold[f, {a}] returns a. Fold[f, {}] remains unevaluated. f may be a
symbol or a pure function; each intermediate application is evaluated
before the next one.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Fold[f, x, {a, b, c, d}]
Out[1]= f[f[f[f[x, a], b], c], d]

In[2]:= Fold[List, x, {a, b, c, d}]
Out[2]= {{{{x, a}, b}, c}, d}

In[3]:= Fold[Times, 1, {a, b, c, d}]
Out[3]= a b c d

In[4]:= Fold[f, {a, b, c, d}]
Out[4]= f[f[f[a, b], c], d]

In[5]:= Fold[{2 #1, 3 #2} &, x, {a, b, c, d}]
Out[5]= {{{{16 x, 24 a}, 12 b}, 6 c}, 3 d}

In[6]:= Fold[f, x, p[a, b, c, d]]
Out[6]= f[f[f[f[x, a], b], c], d]

In[7]:= Fold[2 #1 + #2 &, 0, {1, 0, 1, 1}]
Out[7]= 11

In[8]:= Fold[10 #1 + #2 &, 0, {4, 5, 1, 6, 7, 8}]
Out[8]= 451678
```

## Implementation notes

**Algorithm.** `builtin_fold` (`fold_impl(res, false)`) is a left fold:
`Fold[f, x, {a,b,c}]` → `f[f[f[x,a],b],c]`. The two-argument form `Fold[f,
{a,b,c}]` uses `a` as the seed and folds over the rest; on an empty list this
form stays unevaluated. The seed is pushed into an `ExprBuf` history and the
shared `iter_run` driver is invoked with `fold_step`, which consumes the list
elements in order, computing `apply_binary(f, accumulator, elems[idx++])` (build
`f[acc, e]`, `eval_and_free`). `ebuf_finalize(..., as_list=false)` returns the
last accumulator. `Fold`/`FoldList` are the same `fold_impl`.

**Data structures.** The list's underlying argument array is borrowed (no copy of
the spine); each accumulator value is an owned `Expr*` stored in the history
buffer. The output preserves the input list's head only in the `FoldList` form.

- `Protected`.
- The head of the third argument need not be `List` (any compound expression is accepted).
- `Fold[f, x, {}]` returns `x` (the function is never applied); `Fold[f, {a}]` returns `a`.
- `Fold[f, {}]` remains unevaluated (no seed, no elements).
- Each intermediate application is evaluated before the next one.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Richard Bird, *Introduction to Functional Programming using Haskell*, 2nd ed. (the foldl/reduce accumulator pattern).
- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Fold[f, x, {a, b, c}]
Out[1]= f[f[f[x, a], b], c]
```

```mathematica
In[1]:= Fold[Plus, 0, {1, 2, 3, 4}]
Out[1]= 10
```

```mathematica
In[1]:= Fold[#1*10 + #2 &, 0, {1, 2, 3}]
Out[1]= 123
```

```mathematica
In[1]:= Fold[1/(#2 + #1) &, 0, {1, 1, 1, 1, 1, 1, 1, 1}]
Out[1]= 21/34
```

### Notes

`Fold[f, x, list]` is a left fold: it threads an accumulator (initial value `x`)
through `list`, applying `f[acc, elem]` at each step. The classic `Fold[Plus, 0,
list]` sums a list, while `#1*10 + #2 &` shows the accumulator (`#1`) and current
element (`#2`) being combined to build a number digit by digit. The two-argument
form `Fold[f, list]` uses the list's first element as the seed. `Fold[f, x, {}]`
returns `x` unchanged. The continued-fraction fold `1/(#2 + #1) &` evaluates a
finite continued fraction in exact arithmetic; with eight `1`s it returns the
ratio of consecutive Fibonacci numbers `21/34`, the truncated expansion of the
golden ratio's reciprocal.
