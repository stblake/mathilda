# NestWhile

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
NestWhile[f, expr, test]
    starts with expr and repeatedly applies f while test still yields True.
NestWhile[f, expr, test, m]
    supplies the most recent m results as arguments to test.
NestWhile[f, expr, test, All]
    supplies all results so far as arguments to test.
NestWhile[f, expr, test, {mmin, mmax}]
    delays testing until at least mmin results exist, then passes up to mmax.
NestWhile[f, expr, test, m, max]
    applies f at most max times.
NestWhile[f, expr, test, m, max, n]
    applies f an additional n times after the loop terminates.
NestWhile[f, expr, test, m, max, -n]
    returns the result found when f had been applied n fewer times.

If test[expr] does not yield True initially, NestWhile returns expr.
NestWhile[f, expr, UnsameQ, 2] is equivalent to FixedPoint[f, expr].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= NestWhile[#/2 &, 123456, EvenQ]
Out[1]= 1929

In[2]:= NestWhile[Log, 100., # > 0 &]
Out[2]= -0.859384

In[3]:= NestWhile[Floor[#/2] &, 10, UnsameQ, 2]
Out[3]= 0

In[4]:= NestWhile[#/2 &, 123456, EvenQ, 1, 4]
Out[4]= 7716

In[5]:= NestWhile[Floor[#/2] &, 20, # > 1 &, 1, Infinity]
Out[5]= 1

In[6]:= NestWhile[Floor[#/2] &, 20, # > 1 &, 1, Infinity, 1]
Out[6]= 0

In[7]:= NestWhile[Floor[#/2] &, 20, # > 1 &, 1, Infinity, -1]
Out[7]= 2

In[8]:= NestWhile[# + 1 &, 888, !PrimeQ[#] &]
Out[8]= 907
```

## Implementation notes

- `Protected`.
- If `test[expr]` does not yield `True` initially, the unchanged `expr` is returned.
- Results passed to `test` are in generation order with the most recent last, so e.g. `# > 1 &` inspects the oldest when more than one result is supplied.
- `NestWhile[f, expr, UnsameQ, 2]` is equivalent to `FixedPoint[f, expr]`.
- `NestWhile[f, expr, UnsameQ, All]` continues until any prior value reappears.
- `m` must be a positive integer, `All`, or a 2-element list `{mmin, mmax}` with `1 <= mmin <= mmax` (or `mmax = Infinity`); `max` must be a non-negative integer or `Infinity`; `n` must be an integer. Malformed specs leave `NestWhile` unevaluated.
- Pure functions (`... &`) are supported for both `f` and `test`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
