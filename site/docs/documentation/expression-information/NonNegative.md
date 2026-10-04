# NonNegative

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NonNegative[x]`**

gives True if x is a real number that is positive or zero, and False

<details>
<summary>Notes</summary>

if x is a manifestly negative real number or a non-real complex number. For non-numeric x the expression is left unevaluated. NonNegative is Listable, so it threads over lists element by element.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= NonNegative[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]
Out[1]= {True, True, True, True, False, False, True}

In[2]:= NonNegative[{x, Sin[y]}]
Out[2]= {NonNegative[x], NonNegative[Sin[y]]}

In[3]:= NonNegative[Pi - 3]
Out[3]= True
```

### Applications (5)

Zero is included

```mathematica
In[4]:= NonNegative[0]
Out[4]= True
```

```mathematica
In[5]:= NonNegative[Pi - 3]
Out[5]= True

In[6]:= NonNegative[-3]
Out[6]= False
```

Listable

```mathematica
In[7]:= NonNegative[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]
Out[7]= {True, True, True, True, False, False, True}
```

Non-numeric argument stays symbolic

```mathematica
In[8]:= NonNegative[x]
Out[8]= NonNegative[x]
```

## Implementation notes

**Algorithm.** `builtin_nonnegative` (`src/core.c`) decides whether a numeric
quantity is real and `>= 0`. Like its siblings it first gates on
`is_numeric_quantity` (a non-numeric argument returns `NULL`, leaving the call
unevaluated), then classifies with the shared `numeric_real_sign`: exact
`EXPR_INTEGER` / `EXPR_BIGINT` / `Rational` / `EXPR_REAL` (and `EXPR_MPFR`) are
read directly via `expr_numeric_sign`, everything else is `numericalize`d at
machine precision. A non-real complex value answers `False`; otherwise the verdict
is `sign >= 0`, so zero is included.

**Data structures.** All `Expr`. A packed list or `NDArray` is read in one pass by
`ndint_sign_predicate(res, NDSP_NONNEGATIVE)` (`src/ndinteger.c`) into a `List` of
`True`/`False`; complex buffers, rank > 1, and `Indeterminate` elements are
declined, and `ndarray_delist_and_reeval` falls back to the scalar path reached by
`Listable` threading for an unpacked list.

**Complexity / limits.** `O(1)` per scalar, `O(n)` over a buffer. For inexact or
symbolic-constant inputs the verdict rests on a machine-precision
numericalization. Attributes `Listable`, `Protected`.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [NumericQ](../../expression-information/NumericQ/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_nonnegative.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nonnegative.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`NonNegative[x]` is `True` for a real numeric quantity that is positive **or
zero**; it differs from `Positive` only at zero. A non-real complex value gives
`False`. Exact inputs are decided exactly, inexact and symbolic-constant inputs by
their machine-precision value, and a non-numeric argument is left unevaluated.
`NonNegative` is `Listable` and reads a packed list or `NDArray` straight off the
buffer.
