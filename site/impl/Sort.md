---
source: src/sort.c
---
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
