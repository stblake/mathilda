# NestWhileList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NestWhileList[f, expr, test]
    generates the list {expr, f[expr], f[f[expr]], ...} continuing while
    test applied to the most recent result yields True.
NestWhileList[f, expr, test, m]
    supplies the most recent m results as arguments to test.
NestWhileList[f, expr, test, All]
    supplies all results so far as arguments to test.
NestWhileList[f, expr, test, {mmin, mmax}]
    delays testing until at least mmin results exist, then passes up to mmax.
NestWhileList[f, expr, test, m, max]
    applies f at most max times.
NestWhileList[f, expr, test, m, max, n]
    appends n additional applications of f to the list.
NestWhileList[f, expr, test, m, max, -n]
    drops the last n elements from the list.

NestWhileList[f, expr, UnsameQ, 2] is equivalent to FixedPointList[f, expr].
NestWhileList[f, expr, test, All] is equivalent to
NestWhileList[f, expr, test, {1, Infinity}].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NestWhileList[#/2 &, 123456, EvenQ]
Out[1]= {123456, 61728, 30864, 15432, 7716, 3858, 1929}

In[2]:= NestWhileList[Log, 100., # > 0 &]
Out[2]= {100.0, 4.60517, 1.52718, 0.423423, -0.859384}

In[3]:= NestWhileList[Floor[#/2] &, 20, UnsameQ, 2, 4]
Out[3]= {20, 10, 5, 2, 1}

In[4]:= NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity]
Out[4]= {20, 10, 5, 2, 1}

In[5]:= NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity, 1]
Out[5]= {20, 10, 5, 2, 1, 0}

In[6]:= NestWhileList[Floor[#/2] &, 20, # > 1 &, 1, Infinity, -1]
Out[6]= {20, 10, 5, 2}

In[7]:= NestWhileList[# + 1 &, 899, !PrimeQ[#] &]
Out[7]= {899, 900, 901, 902, 903, 904, 905, 906, 907}

In[8]:= NestWhileList[Mod[2 #, 19] &, 2, # != 1 &]
Out[8]= {2, 4, 8, 16, 13, 7, 14, 9, 18, 17, 15, 11, 3, 6, 12, 5, 10, 1}
```

## Implementation notes

- `Protected`.
- Results are listed in generation order, including the final element on which `test` yielded a non-`True` value (or the last element produced when `max` iterations were reached).
- If `test[expr]` does not yield `True` initially, the result is just `{expr}`.
- `NestWhileList[f, expr, UnsameQ, 2]` is equivalent to `FixedPointList[f, expr]`.
- `NestWhileList[f, expr, test, All]` is equivalent to `NestWhileList[f, expr, test, {1, Infinity}]`.
- `NestWhileList[f, expr, UnsameQ, All]` continues until a previously-seen value reappears, and the repeat is included as the last element of the list.
- `m` must be a positive integer, `All`, or a 2-element list `{mmin, mmax}` with `1 <= mmin <= mmax` (or `mmax = Infinity`); `max` must be a non-negative integer or `Infinity`; `n` must be an integer. Malformed specs leave `NestWhileList` unevaluated.
- Pure functions (`... &`) are supported for both `f` and `test`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
