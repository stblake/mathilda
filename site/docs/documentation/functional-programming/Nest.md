# Nest

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Nest[f, expr, n]
    gives an expression with f applied n times to expr.

n must be a non-negative integer. Nest[f, expr, 0] returns expr. The
function f may be a symbol or a pure function. Each iteration evaluates
f applied to the current value before proceeding.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Nest[f, x, 3]
Out[1]= f[f[f[x]]]

In[2]:= Nest[(1 + #)^2 &, 1, 3]
Out[2]= 676

In[3]:= Nest[(1 + #)^2 &, x, 5]
Out[3]= (1 + (1 + (1 + (1 + (1 + x)^2)^2)^2)^2)^2

In[4]:= Nest[Sqrt, 100.0, 4]
Out[4]= 1.33352

In[5]:= Nest[1/(1 + #) &, x, 5]
Out[5]= 1/(1 + 1/(1 + 1/(1 + 1/(1 + 1/(1 + x)))))

In[6]:= Nest[x^# &, x, 6]
Out[6]= x^x^x^x^x^x^x

In[7]:= Nest[#(1 + 0.05) &, 1000, 10]
Out[7]= 1628.89

In[8]:= Nest[(# + 2/#)/2 &, 1.0, 5]
Out[8]= 1.41421
```

## Implementation notes

`builtin_nest` (via `nest_impl(res, false)`) applies `f` to `expr` exactly `n`
times and returns the final value; `n` must be a non-negative integer or the
call stays unevaluated. It seeds a growable history buffer `ExprBuf` with a copy
of `expr` and drives the shared generic runner `iter_run` with `nest_step`, which
on each step computes `apply_unary(f, last)` (build `f[last]`, `eval_and_free`).
`ebuf_finalize(..., as_list=false)` returns the last history element and frees the
rest. `Nest` and `NestList` are the same `nest_impl`, differing only in the
`as_list` flag.

- `Protected`.
- `n` must be a non-negative integer; `Nest[f, expr, 0]` returns `expr` unchanged.
- The function `f` may be a symbol, a built-in, or a pure function (`... &`).
- Each iteration evaluates `f[current]` before proceeding, so numeric computations collapse immediately.
- Returns unevaluated if `n` is not a non-negative integer or the argument count is wrong.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Nest[f, x, 3]
Out[1]= f[f[f[x]]]
```

```mathematica
In[1]:= Nest[#^2 &, 2, 3]
Out[1]= 256
```

```mathematica
In[1]:= Nest[1/(1 + #) &, x, 2]
Out[1]= 1/(1 + 1/(1 + x))
```

### Notes

`Nest[f, expr, n]` applies `f` to `expr` exactly `n` times and returns only the
final result. `n` must be a non-negative integer; `Nest[f, expr, 0]` returns
`expr` unchanged. Each intermediate application is evaluated before the next, so
numeric iterations like `#^2 &` collapse to a single number (`2 -> 4 -> 16 ->
256`). Use `NestList` instead when the intermediate values are also wanted.
