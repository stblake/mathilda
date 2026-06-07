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

**Algorithm.** `builtin_nestwhile` (`nestwhile_impl(res, false)`) iterates `f`
starting from `expr` while a predicate `test` holds, returning the last value. It
accepts up to six arguments `(f, expr, test, m, max, n)`. `m` controls how many
of the most-recent history entries are passed to `test` (an integer, `All`, or
`{mmin, mmax|Infinity}`; default `1`); `max` caps the number of `f`-applications
(integer or `Infinity`); `n` is post-processing — positive `n` applies `f` that
many extra times, negative `n` drops `|n|` iterates from the end.

The shared `iter_run` driver runs `nestwhile_step`: once at least `m_min` history
entries exist it builds `test[recent...]` from the last `min(count, m_max)`
entries, evaluates it, and halts (without applying `f` again) if the result is
not `True`; otherwise it appends `apply_unary(f, last)`. For an unbounded `max`,
`iter_run` is given `ITER_SAFETY_CAP` (1,000,000) as a safety limit, returning
`NULL` if exceeded. Malformed argument specs return `NULL` (unevaluated).
Post-processing reuses `nest_step` (positive `n`) or `ebuf_truncate` (negative
`n`).

**Data structures.** Growable `ExprBuf` history of owned `Expr*`; the trailing
window passed to `test` is a freshly copied argument array per step.

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

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
