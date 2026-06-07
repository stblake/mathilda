# Sort

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Sort[list] sorts the elements of list into canonical order.
Sort[list, p] sorts using the ordering function p.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Sort[{d, b, c, a}]
Out[1]= {a, b, c, d}

In[2]:= Sort[{Pi, E, 2, 3, 1, Sqrt[2]}, Less]
Out[2]= {1, Sqrt[2], 2, E, 3, Pi}
```

## Implementation notes

**Algorithm.** `builtin_sort` deep-copies the argument list's elements into an `Expr**` array
and sorts it in place with the C library `qsort`. With no ordering function it uses
`expr_compare` (the canonical Mathematica `Order`); with a second argument `p` it calls
`p[a, b]`, evaluates the result, and treats `True`/`1` as "in order" and `False`/`-1` as "out
of order". The custom comparator is passed to `qsort` through a file-scope `current_sort_p`
pointer (saved/restored around the call so reentrant sorts nest correctly). The original head is
preserved, so `Sort` works on any expression, not just `List`.

**Canonical order (`expr_compare`).** Defined in this file (co-located with `Sort`/`OrderedQ`
and also used by the evaluator's `Orderless` argument-sorting): (1) numeric atoms (Integer,
Real, Rational, BigInt, MPFR) sort first by value, integers compared exactly via GMP; (2)
strings next, case-insensitive then case-sensitive lexicographic; (3) everything else is
compared by a **polynomial degree vector** — collect every symbol name in either operand, sort
those names in reverse-alphabetical order (so the lex-last variable is most significant), and
lexicographically compare the per-variable degrees (`expr_poly_degree`, which returns +∞ for
non-polynomial occurrences). This gives the grevlex-with-reverse-alpha display order. (4) Ties
break structurally: bare symbol before compound, then head, arity, and args recursively.

**Complexity.** `O(n log n)` comparisons; each `expr_compare` is itself `O(symbols × tree
size)` because it rebuilds the symbol set per pair. The comparator is deliberately made
symmetric (order-independent symbol collection) so `qsort` cannot oscillate on `Orderless`
heads with many unknowns.

- `Protected`.
- Uses an efficient quicksort algorithm.
- Canonical order:
    - Real numbers by numerical value.
    - Complex numbers by real part, then imaginary part magnitude.
    - Strings in dictionary order (lowercase before uppercase).
    - Symbols by name.
    - Expressions by length, then head, then parts depth-first.
- Polynomial order: `x^n` sorts relative to `x`.
- Numeric coefficient stripping: when a `Times` term has a leading numeric factor (Integer, Real, BigInt, MPFR, `Rational[n,d]`, `Complex[re,im]`, or a radical such as `Sqrt[2] = Power[2, 1/2]`), that factor is ignored when computing the term's main factor. As a result, `1 + x^2 + Sqrt[2] x` is canonicalised to `1 + Sqrt[2] x + x^2`, matching Mathematica's main-factor-first ordering.
- Custom ordering function `p` can return `1`, `0`, `-1`, `True`, or `False`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
