# NestList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NestList[f, expr, n]
    gives a list of the results of applying f to expr 0 through n times.

The result is a list of length n+1 whose first element is expr and
whose (k+1)-th element is f applied k times to expr. n must be a
non-negative integer. f may be a symbol or a pure function; each
intermediate application is evaluated before the next one.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NestList[f, x, 4]
Out[1]= {x, f[x], f[f[x]], f[f[f[x]]], f[f[f[f[x]]]]}

In[2]:= NestList[Cos, 1.0, 10]
Out[2]= {1.0, 0.540302, 0.857553, 0.65429, 0.79348, 0.701369, 0.76396, 0.722102, 0.750418, 0.731404, 0.744237}

In[3]:= NestList[(1 + #)^2 &, x, 3]
Out[3]= {x, (1 + x)^2, (1 + (1 + x)^2)^2, (1 + (1 + (1 + x)^2)^2)^2}

In[4]:= NestList[Sqrt, 100.0, 4]
Out[4]= {100.0, 10.0, 3.16228, 1.77828, 1.33352}

In[5]:= NestList[2 # &, 1, 10]
Out[5]= {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024}

In[6]:= NestList[# + 1 &, 0, 10]
Out[6]= {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

In[7]:= NestList[#^2 &, 2, 6]
Out[7]= {2, 4, 16, 256, 65536, 4294967296, 18446744073709551616}

In[8]:= NestList[#(1 + 0.05) &, 1000, 10]
Out[8]= {1000, 1050.0, 1102.5, 1157.62, 1215.51, 1276.28, 1340.1, 1407.1, 1477.46, 1551.33, 1628.89}
```

## Implementation notes

`builtin_nestlist` is `nest_impl(res, true)`: it returns the full iteration
history `{expr, f[expr], ..., Nest[f,expr,n]}`. It seeds an `ExprBuf` with a copy
of `expr` and runs the shared `iter_run` driver with `nest_step` (each step is
`apply_unary(f, last)` = build `f[last]` and `eval_and_free`). With
`as_list=true`, `ebuf_finalize` wraps the entire kept history in a `List` head.
Shares all machinery with `Nest`; `n` must be a non-negative integer.

- `Protected`.
- Returns a list of length `n + 1` whose first element is `expr` and whose `(k+1)`-th element is `f` applied `k` times to `expr`.
- `n` must be a non-negative integer; `NestList[f, expr, 0]` returns `{expr}`.
- The function `f` may be a symbol, a built-in, or a pure function (`... &`).
- Each iteration evaluates `f[current]` before proceeding, so numeric computations collapse immediately.
- Returns unevaluated if `n` is not a non-negative integer or the argument count is wrong.
- `Last[NestList[f, expr, n]]` is equivalent to `Nest[f, expr, n]`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= NestList[f, x, 3]
Out[1]= {x, f[x], f[f[x]], f[f[f[x]]]}
```

```mathematica
In[1]:= NestList[2 # &, 1, 5]
Out[1]= {1, 2, 4, 8, 16, 32}
```

The successive convergents of Newton's iteration for `Sqrt[2]`, kept exact —
each entry roughly doubles the number of correct digits:

```mathematica
In[1]:= NestList[(# + 2/#)/2 &, 1, 4]
Out[1]= {1, 3/2, 17/12, 577/408, 665857/470832}
```

Building the first few convergents of a continued fraction symbolically:

```mathematica
In[1]:= NestList[1/(1 + #) &, x, 3]
Out[1]= {x, 1/(1 + x), 1/(1 + 1/(1 + x)), 1/(1 + 1/(1 + 1/(1 + x)))}
```

### Notes

`NestList[f, expr, n]` returns a list of length `n + 1`: the first element is
`expr` and each subsequent element applies one more `f`. It is the
value-collecting companion of `Nest`, useful for capturing an entire orbit or
sequence (powers of two, a counter, the iterates of a map). `n` must be a
non-negative integer, and each intermediate value is evaluated before being
appended.
