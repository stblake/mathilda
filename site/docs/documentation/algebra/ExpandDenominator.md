# ExpandDenominator

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ExpandDenominator[expr]`**

expands out products and powers that appear as denominators in expr.

<details>
<summary>Notes</summary>

ExpandDenominator works only on negative integer powers. ExpandDenominator applies only to the top level in expr. ExpandDenominator leaves the numerator unexpanded. ExpandDenominator automatically threads over lists, as well as equations, inequalities, and logic functions.

</details>

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ExpandDenominator[(x-1)(x-2)/((x-3)(x-4))]
Out[1]= ((-2 + x) (-1 + x))/(12 - 7 x + x^2)

In[2]:= ExpandDenominator[1/(x+1) + 2/(x+1)^2 + 3/(x+1)^3]
Out[2]= 1/(1 + x) + 2/(1 + 2 x + x^2) + 3/(1 + 3 x + 3 x^2 + x^3)

In[3]:= ExpandDenominator[(a+b)(a-b)/((c+d)(c-d))]
Out[3]= ((a + b) (a - b))/(c^2 - d^2)
```

### Applications (3)

```mathematica
In[4]:= ExpandDenominator[(a + b)^2 / (c + d)^2]
Out[4]= (a + b)^2/(c^2 + 2 c d + d^2)

In[5]:= ExpandDenominator[1 / ((x + 1)(x + 2)(x + 3))]
Out[5]= 1/(6 + 11 x + 6 x^2 + x^3)

In[6]:= ExpandDenominator[((x + 1)(x + 2)) / (y (y + 1))]
Out[6]= ((1 + x) (2 + x))/(y + y^2)
```

## Implementation notes

**Algorithm.** `builtin_expand_denominator` (in `src/expand.c`) calls `expr_expand_denominator`, the mirror of `ExpandNumerator`: it expands only the denominator. For a `Times` it collects the negative-integer-power factors (`is_negative_int_power`), rebuilds them as positive powers, multiplies them into a single denominator product, runs `expr_expand` on that product, and re-inverts it as `Power[expandedDen, -1]` multiplied back against the untouched numerator factors. A bare `Power[base, -k]` is handled by expanding `base^k` and re-inverting. Threading over `List`, equations, inequalities, logic heads and `Plus` matches `ExpandNumerator`.

**Data structures.** Parallel `Expr**` buffers for numerator factors and positive-power denominator factors; only the denominator product is expanded, the numerator is copied through.

- `Protected`.
- Acts only on factors with negative integer exponents (the "denominator part" of `expr`).
- Combines all denominator factors of a top-level `Times` into a single expanded polynomial wrapped in `Power[..., -1]`.
- Applies only to the top level in `expr`; it does not descend into function bodies.
- Leaves the numerator factors unchanged.
- Threads over `List`, `Equal`, `Unequal`, `Less`, `LessEqual`, `Greater`, `GreaterEqual`, `And`, `Or`, `Not`, and `Plus`.

**Attributes:** `Protected`.

## References

**See also:** [Times](../../arithmetic/Times/), [List](../../other-advanced/List/), [Equal](../../comparisons/Equal/), [Unequal](../../comparisons/Unequal/), [Less](../../comparisons/Less/), [LessEqual](../../comparisons/LessEqual/), [Greater](../../comparisons/Greater/), [GreaterEqual](../../comparisons/GreaterEqual/)

- Source: [`src/expand.c`](https://github.com/stblake/mathilda/blob/main/src/expand.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_expandfrac.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expandfrac.c)

## Notes & additional examples

### Notes

`ExpandDenominator` multiplies out only the denominator (the negative integer
powers), leaving the numerator factored — the mirror image of
`ExpandNumerator`. The product of three linear factors above expands to the
cubic `6 + 11 x + 6 x^2 + x^3` in the denominator. It threads over lists,
equations, inequalities, and logic functions.
