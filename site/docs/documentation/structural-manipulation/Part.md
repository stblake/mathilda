# Part

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

expr\[\[i\]\] or Part\[expr, i\] gives the i-th part of expr. expr\[\[-i\]\] counts from the end. expr\[\[0\]\] gives the head of expr. expr\[\[i, j, ...\]\] or Part\[expr, i, j, ...\] is equivalent to expr\[\[i\]\]\[\[j\]\]..., descending into nested parts. expr\[\[{i1, i2, ...}\]\] gives a list of the parts i1, i2, ... of expr (wrapped in the head of expr). expr\[\[m;;n\]\] / expr\[\[m;;n;;s\]\] gives the span of parts m through n (with optional step s); ;; alone or All means all parts. Part is treated as atomic on Integer, Real, String, Symbol, Rational\[n, d\], and Complex\[re, im\]; Part\[atom, i\] for i != 0 stays unevaluated. Indices are 1-based and may be negative; out-of-range indices leave the expression unevaluated.

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= {a, b, c, d}[[2]]
Out[1]= b

In[2]:= {a, b, c, d}[[-1]]
Out[2]= d

In[3]:= 123[[0]]
Out[3]= Integer
```

### Applications (10)

Extract a single part (1-based), count from the end with a negative index, and
get the head with index `0`:

```mathematica
In[1]:= {a, b, c, d}[[2]]
Out[1]= b

In[2]:= {a, b, c, d}[[-1]]
Out[2]= d

In[3]:= {a, b, c, d}[[0]]
Out[3]= List
```

Multi-index parts descend into nested structure, and a list of indices gathers
several parts at once:

```mathematica
In[1]:= {{1, 2}, {3, 4}}[[2, 1]]
Out[1]= 3

In[2]:= {a, b, c, d}[[{1, 3}]]
Out[2]= {a, c}
```

Spans `m ;; n ;; s` slice with an optional step:

```mathematica
In[1]:= {a, b, c, d, e, f}[[1 ;; 6 ;; 2]]
Out[1]= {a, c, e}
```

`Part` is the workhorse of matrix manipulation. Combine `All` with an index to
pull out a whole column, and a pair of index lists to carve out an arbitrary
submatrix — here the corner 2×2 block of a 3×3 matrix, whose trace is the sum of
the corner entries:

```mathematica
In[1]:= m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

In[2]:= m[[All, 2]]
Out[2]= {2, 5, 8}

In[3]:= Tr[m[[{1, 3}, {1, 3}]]]
Out[3]= 10
```

Because everything is an expression, `Part` reaches into non-list heads too —
the second summand of a sum:

```mathematica
In[1]:= (a + b + c)[[2]]
Out[1]= b
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Sort 4x10^6 | 42.2 s | 68.7 s | 111 s |
| gather v[[idx]], 4x10^6 | 16.8 s | 6.66 s | 7.18 s |
| Union of 4x10^6 integers | 12.4 s | 71.1 s | 376 s |
| Reverse 4x10^6 | 5.37 s | 0.297 s | 0.982 s |
| Join two 2x10^6 | 0.899 s | 0.6 s | 0.397 s |
| RotateLeft 4x10^6 by 1000 | 0.855 s | 0.307 s | 0.456 s |

## Implementation notes

**Algorithm.** `builtin_part` extracts elements by index path, delegating to the recursive
`expr_part(expr, indices, nindices)`. Each index level may be: a positive or negative integer
(`-k` resolves to `len + k + 1`); `0`, which extracts the head and is allowed even on atoms; a
`List` of indices (extract several, returning a list); `All`; or a `Span` (`i;;j;;k`) built by
the parser from `;;` syntax, which `expr_part` resolves into an explicit element range with the
given start/end/step (negative endpoints wrap, `UpTo`/`All` endpoints clamp). Index paths apply
left to right, descending one structural level per index. Out-of-range or non-integer indices on
atoms yield `NULL` (unevaluated).

- Supports negative indices to count from the end (`-1` is the last element).
- `expr[[0]]` returns the `Head` of the expression. This is permitted even for atomic expressions.
- Mapping `All` across a dimension allows column extraction from matrices.

**Attributes:** `NHoldRest`, `Protected`.

## See also

[Head](../../structural-manipulation/Head/)

## References

- Source: [`src/part.c`](https://github.com/stblake/mathilda/blob/main/src/part.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_matsol_methods.c`](https://github.com/stblake/mathilda/blob/main/tests/test_matsol_methods.c)
- Tests: [`tests/test_nresidue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nresidue.c)
- Tests: [`tests/test_nroots.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nroots.c)

## Notes & additional examples

### Notes

`expr[[i]]` (= `Part[expr, i]`) is 1-based; `expr[[-i]]` counts from the end and
`expr[[0]]` returns the head. Multi-index `expr[[i, j, ...]]` descends through
nested parts, a list of indices `expr[[{i1, i2, ...}]]` gathers several (rewrapped
in the original head), and spans `m ;; n ;; s` (with `All` or `;;` meaning "every
part") slice ranges. Part operates on any head, not just `List`, so it indexes
sums, products, and arbitrary symbolic structure uniformly. It is atomic on
Integer, Real, String, Symbol, Rational, and Complex; out-of-range or
unsupported indices leave the expression unevaluated rather than erroring.
