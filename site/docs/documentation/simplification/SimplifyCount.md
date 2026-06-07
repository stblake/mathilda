# SimplifyCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SimplifyCount[expr]
    The complexity measure used by Simplify when no
    ComplexityFunction option (or ComplexityFunction -> Automatic) is
    given. Counts subexpressions; integers contribute their decimal
    digit count plus a constant for the sign. Real numbers contribute
    2 (NumberQ but not Integer/Rational).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= SimplifyCount[100 Log[2]]
Out[1]= 6

In[2]:= SimplifyCount[Log[2^100]]
Out[2]= 32

In[3]:= SimplifyCount[1/2]
Out[3]= 3

In[4]:= SimplifyCount[3.14]
Out[4]= 2
```

## Implementation notes

**Algorithm.** `builtin_simplify_count` returns `simp_default_complexity(arg)` as
an Integer — the complexity measure Simplify uses by default (when no
`ComplexityFunction`, or `ComplexityFunction -> Automatic`, is given). It is a
recursive leaf count with Mathematica-faithful adjustments: `Symbol -> 1`;
`String -> 1`; `Real`/`MPFR -> 2` (NumberQ but not Integer/Rational); `Integer 0
-> 1`, positive integer `p -> digits(p)`, negative `-> digits(|p|) + 1` (the
leading minus sign costs one unit, computed by `int_digit_count_int64` /
`mpz_sizeinbase` for bigints); `Rational[n,d]` and `Complex[re,im]` add `1` for
the wrapper plus the counts of their two parts; any other `Function` is
`count(head) + sum count(args)`. The negative-integer and digit-count rules are
what keep e.g. `100 Log[2]` (score 6) preferred over `Log[2^100]` (score 32),
preventing Simplify from folding into giant exact integers.

**Data structures.** Pure recursive walk over the `Expr*` tree; no allocation.
The same `size_t simp_default_complexity` function is consumed internally by
`score_with_func`. `ATTR_LISTABLE`, so `SimplifyCount[{a, b}]` threads.

- `Listable`, `Protected`.
- Per node: a symbol, the integer `0`, or a string counts `1`; a positive integer

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/simp/simp_complexity.c`](https://github.com/stblake/mathilda/blob/main/src/simp/simp_complexity.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)
