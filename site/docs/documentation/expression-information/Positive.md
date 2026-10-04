# Positive

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Positive[x]`**

gives True if x is a positive real number, and False if x is a

<details>
<summary>Notes</summary>

manifestly negative real number, a non-real complex number, or zero. For non-numeric x the expression is left unevaluated. Positive is Listable, so it threads over lists element by element.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Positive[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]
Out[1]= {True, True, True, False, False, False, True}

In[2]:= Positive[{x, Sin[y]}]
Out[2]= {Positive[x], Positive[Sin[y]]}

In[3]:= Positive[Sqrt[-2]]
Out[3]= False
```

### Applications (5)

```mathematica
In[4]:= Positive[3]
Out[4]= True
```

A symbolic constant is decided by its numeric value

```mathematica
In[5]:= Positive[Pi - 3]
Out[5]= True
```

A non-real complex value is not positive

```mathematica
In[6]:= Positive[Sqrt[-2]]
Out[6]= False
```

Listable: threads over the list

```mathematica
In[7]:= Positive[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]
Out[7]= {True, True, True, False, False, False, True}
```

Non-numeric argument stays symbolic

```mathematica
In[8]:= Positive[x]
Out[8]= Positive[x]
```

## Implementation notes

**Algorithm.** `builtin_positive` (`src/core.c`) decides the sign of a numeric
quantity. It first gates on `is_numeric_quantity`: a non-numeric argument (a bare
symbol, `Positive[x]`) returns `NULL`, so the call is left unevaluated and the
symbolic expression flows on. For a numeric argument the shared helper
`numeric_real_sign` classifies reality and sign — exact `EXPR_INTEGER` /
`EXPR_BIGINT` / `Rational` / `EXPR_REAL` (and `EXPR_MPFR`) are read directly with
`expr_numeric_sign`, everything else is `numericalize`d at machine precision
(`numeric_machine_spec()`). A value that numericalizes to a genuinely non-real
complex number answers `False`; otherwise the verdict is `sign > 0`.

**Data structures.** Everything is `Expr`. A packed list or `NDArray` takes a
one-pass fast path: `ndint_sign_predicate(res, NDSP_POSITIVE)` (`src/ndinteger.c`)
reads the machine buffer and emits a `List` of `True`/`False` with no per-element
`Expr` and no evaluator round-trip. It declines the shapes it cannot cover —
complex buffers, rank > 1, and `Indeterminate` elements (unordered under IEEE) —
and `ndarray_delist_and_reeval` then threads the ordinary scalar path, which is
reached anyway by the `Listable` attribute for an unpacked list.

**Complexity / limits.** `O(1)` per scalar, `O(n)` element-wise over a buffer.
The verdict for an inexact or symbolic-constant argument rests on a
machine-precision numericalization, so a quantity indistinguishable from zero at
machine precision is classified by its machine sign. Attributes `Listable`,
`Protected`.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [NumericQ](../../expression-information/NumericQ/), [Negative](../../expression-information/Negative/), [NonNegative](../../expression-information/NonNegative/), [NonPositive](../../expression-information/NonPositive/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_positive.c`](https://github.com/stblake/mathilda/blob/main/tests/test_positive.c)

## Notes & additional examples

### Notes

`Positive[x]` is `True` only for a real, strictly positive numeric quantity.
Zero gives `False`, and so does any non-real complex value — `Positive` is a test
of a real sign, not of "has a positive real part". Exact integers, rationals, and
bigints are decided exactly; reals, symbolic constants, and numeric-function calls
are classified by their machine-precision numeric value.

A non-numeric argument (one for which `NumericQ` is `False`) is left unevaluated,
so symbolic expressions flow through the evaluator unchanged. `Positive` is
`Listable`, and a packed list or `NDArray` is read straight off the buffer; the
result is a list of `True`/`False`, which no buffer holds, so it comes back as an
ordinary list.
