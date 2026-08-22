# ExpandNumerator

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ExpandNumerator[expr]`**

expands out products and powers that appear in the numerator of expr.

<details>
<summary>Notes</summary>

ExpandNumerator works on terms that have positive integer exponents. ExpandNumerator applies only to the top level in expr. ExpandNumerator does not separate the fraction; Expand does. ExpandNumerator leaves the denominator unexpanded. ExpandNumerator automatically threads over lists, as well as equations, inequalities, and logic functions.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ExpandNumerator[(x-1)(x-2)/((x-3)(x-4))]
Out[1]= (2 - 3 x + x^2)/((-4 + x) (-3 + x))

In[2]:= ExpandNumerator[(a+b)^2/x + (c+d)(c-d)/y]
Out[2]= (a^2 + 2 a b + b^2)/x + (c^2 - d^2)/y

In[3]:= ExpandNumerator[x == (a+b)^2/c && y >= (a-b)^2/c]
Out[3]= x == (a^2 + 2 a b + b^2)/c && y >= (a^2 - 2 a b + b^2)/c
```

### Applications (3)

```mathematica
In[4]:= ExpandNumerator[(a + b)^2 / (c + d)^2]
Out[4]= (a^2 + 2 a b + b^2)/(c + d)^2

In[5]:= ExpandNumerator[((x + 1)(x + 2)) / (y (y + 1))]
Out[5]= (2 + 3 x + x^2)/(y (1 + y))

In[6]:= ExpandNumerator[(1 + x)^3 / x^2 == (1 + y)^2 / y]
Out[6]= (1 + 3 x + 3 x^2 + x^3)/x^2 == (1 + 2 y + y^2)/y
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

## References

**See also:** [Expand](../../algebra/Expand/), [List](../../other-advanced/List/), [Equal](../../comparisons/Equal/), [Unequal](../../comparisons/Unequal/), [Less](../../comparisons/Less/), [LessEqual](../../comparisons/LessEqual/), [Greater](../../comparisons/Greater/), [GreaterEqual](../../comparisons/GreaterEqual/)

- Source: [`src/expand.c`](https://github.com/stblake/mathilda/blob/main/src/expand.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_expandfrac.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expandfrac.c)

## Notes & additional examples

### Notes

`ExpandNumerator` distributes products and positive integer powers in the
numerator only, leaving the denominator factored. Unlike `Expand`, it does not
split the fraction into separate terms. It threads over lists and relations, so
both sides of the equation above have their numerators expanded independently.
