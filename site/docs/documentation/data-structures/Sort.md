# Sort

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Sort[list] sorts the elements of list into canonical order.`**

**`Sort[list, p] sorts using the ordering function p.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= Sort[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[1]= <|"b" -> 1, "c" -> 2, "a" -> 3|>

In[2]:= SortBy[<|"a" -> {9}, "b" -> {1}|>, First]
Out[2]= <|"b" -> {1}, "a" -> {9}|>

In[3]:= Total[<|"a" -> 3, "b" -> 1, "c" -> 2|>]
Out[3]= 6

In[4]:= Join[<|"a" -> 1, "b" -> 2|>, <|"b" -> 3, "c" -> 4|>]
Out[4]= <|"a" -> 1, "b" -> 3, "c" -> 4|>
```

### Applications (6)

```mathematica
In[5]:= Sort[{3, 1, 2}]
Out[5]= {1, 2, 3}

In[6]:= Sort[{5, 3, 8, 1}, Greater]
Out[6]= {8, 5, 3, 1}

In[7]:= Sort[{x^2, x, 1, x^3}]
Out[7]= {1, x, x^2, x^3}

In[8]:= Sort[{"banana", "apple", "cherry"}]
Out[8]= {"apple", "banana", "cherry"}

In[9]:= Sort[{{2, 1}, {1, 3}, {1, 2}}]
Out[9]= {{1, 2}, {1, 3}, {2, 1}}

In[10]:= Sort[Range[10], (Mod[#1, 3] < Mod[#2, 3]) &]
Out[10]= {3, 9, 6, 10, 1, 7, 4, 2, 8, 5}
```

## Performance

Measured on arm64 Darwin at commit `2dea9cc05`.

| case | n | time |
|---|---:|---:|
| list of machine reals | 1,000 | 37 us |
| list of machine reals | 10,000 | 117 us |
| list of machine reals | 100,000 | 1.2 ms |

## Implementation notes

**Algorithm.** `builtin_sort` deep-copies the argument list's elements into an `Expr**` array
and sorts it in place with the C library `qsort`. With no ordering function it uses
`expr_compare` (the canonical `Order`); with a second argument `p` it calls
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

**Attributes:** `Protected`.

## References

**See also:** [SortBy](../../data-structures/SortBy/), [Total](../../arithmetic/Total/), [Min](../../data-structures/Min/), [Max](../../data-structures/Max/), [Join](../../data-structures/Join/)

- Source: [`src/sort.c`](https://github.com/stblake/mathilda/blob/main/src/sort.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_bignum_rational_numeric.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bignum_rational_numeric.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)

## Notes & additional examples

### Notes

`Sort[list]` orders elements by Mathilda's canonical ordering, which compares numbers numerically, strings lexicographically, and structured expressions component-by-component (so the nested lists sort by first element, then second). `Sort[list, p]` uses an ordering predicate `p[a, b]` instead: `Greater` reverses to descending order, and a pure function such as `Mod[#1, 3] < Mod[#2, 3]` groups by residue class. The sort is stable, so equal-ranked elements keep their original relative order.
