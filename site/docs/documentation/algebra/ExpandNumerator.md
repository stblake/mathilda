# ExpandNumerator

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ExpandNumerator[expr]
    expands out products and powers that appear in the numerator of expr.
ExpandNumerator works on terms that have positive integer exponents.
ExpandNumerator applies only to the top level in expr.
ExpandNumerator does not separate the fraction; Expand does.
ExpandNumerator leaves the denominator unexpanded.
ExpandNumerator automatically threads over lists, as well as equations,
    inequalities, and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= ExpandNumerator[(x-1)(x-2)/((x-3)(x-4))]
Out[1]= (2 - 3 x + x^2)/((-4 + x) (-3 + x))

In[2]:= ExpandNumerator[(a+b)^2/x + (c+d)(c-d)/y]
Out[2]= (a^2 + 2 a b + b^2)/x + (c^2 - d^2)/y

In[3]:= ExpandNumerator[x == (a+b)^2/c && y >= (a-b)^2/c]
Out[3]= x == (a^2 + 2 a b + b^2)/c && y >= (a^2 - 2 a b + b^2)/c
```

## Implementation notes

**Algorithm.** `builtin_expand_numerator` (in `src/expand.c`) calls `expr_expand_numerator`, which separates an expression's numerator from its denominator and expands only the numerator. For a `Times`, it partitions factors into denominator factors (those of the form `Power[base, k]` with `k` a negative integer, detected by `is_negative_int_power`) and the rest; the non-denominator product is run through `expr_expand`, then recombined with the untouched denominator factors. A bare `Power` with negative integer exponent is a pure denominator and is returned unchanged; a positive/symbolic power is expanded at the top level. It threads over `List`, equations, inequalities, `And`/`Or`/`Not`, and `Plus` (the `is_thread_head` set), expanding per-summand.

**Data structures.** Separate `Expr**` accumulators for numerator and denominator factors; the denominator is preserved verbatim while only the numerator passes through `expr_expand`.

- `Protected`.
- Acts only on factors with positive integer exponents (the "numerator part" of `expr`).
- Applies only to the top level in `expr`; it does not descend into function bodies.
- Leaves the denominator factors (those with negative integer exponents) unchanged.
- Does not separate the fraction into a sum of fractions; only `Expand` does that.
- Threads over `List`, `Equal`, `Unequal`, `Less`, `LessEqual`, `Greater`, `GreaterEqual`, `And`, `Or`, `Not`, and `Plus` (so each summand of a sum-of-fractions is processed independently).

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/expand.c`](https://github.com/stblake/mathilda/blob/main/src/expand.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
