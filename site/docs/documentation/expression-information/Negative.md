# Negative

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Negative[x]`**

gives True if x is a negative real number, and False if x is a

<details>
<summary>Notes</summary>

manifestly non-negative real number (including zero) or a non-real complex number. For non-numeric x the expression is left unevaluated. Negative is Listable, so it threads over lists element by element.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Negative[{1.6, 3/4, Pi, 0, -5, 1 + I, Sin[10^5]}]
Out[1]= {False, False, False, False, True, False, False}

In[2]:= Negative[{x, Sin[y]}]
Out[2]= {Negative[x], Negative[Sin[y]]}

In[3]:= Negative[1 - Pi]
Out[3]= True
```

### Applications (5)

```mathematica
In[4]:= Negative[1 - Pi]
Out[4]= True

In[5]:= Negative[-5]
Out[5]= True
```

A non-real complex value is not negative

```mathematica
In[6]:= Negative[1 + I]
Out[6]= False
```

Listable: threads over the list

```mathematica
In[7]:= Negative[{1.6, 3/4, Pi, 0, -5}]
Out[7]= {False, False, False, False, True}
```

Non-numeric argument stays symbolic

```mathematica
In[8]:= Negative[x]
Out[8]= Negative[x]
```

## Implementation notes

**Algorithm.** `builtin_negative` (`src/core.c`) is the sign-mirror of
`builtin_positive`. It gates on `is_numeric_quantity` — a non-numeric argument
returns `NULL` and stays unevaluated — then calls the shared `numeric_real_sign`
helper, which reads exact `EXPR_INTEGER` / `EXPR_BIGINT` / `Rational` /
`EXPR_REAL` (and `EXPR_MPFR`) directly via `expr_numeric_sign` and
`numericalize`s anything else at machine precision. A value that numericalizes to
a non-real complex number answers `False`; otherwise the verdict is `sign < 0`.

**Data structures.** All `Expr`. A packed list or `NDArray` is read in one pass by
`ndint_sign_predicate(res, NDSP_NEGATIVE)` (`src/ndinteger.c`), emitting a `List`
of `True`/`False` straight off the buffer; it declines complex buffers, rank > 1,
and `Indeterminate` elements, after which `ndarray_delist_and_reeval` falls back
to the scalar path that `Listable` threading reaches for an unpacked list.

**Complexity / limits.** `O(1)` per scalar, `O(n)` over a buffer. For inexact or
symbolic-constant inputs the decision rests on a machine-precision
numericalization, so a value indistinguishable from zero at machine precision is
classified by its machine sign. Attributes `Listable`, `Protected`.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [NumericQ](../../expression-information/NumericQ/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_negative.c`](https://github.com/stblake/mathilda/blob/main/tests/test_negative.c)

## Notes & additional examples

### Notes

`Negative[x]` is `True` only for a real, strictly negative numeric quantity. Zero
gives `False` (use `NonPositive` to include it), and so does any non-real complex
value. Exact integers, rationals, and bigints are decided exactly; everything else
numeric is classified by its machine-precision value.

A non-numeric argument is left unevaluated. `Negative` is `Listable` and reads a
packed list or `NDArray` straight off the buffer, returning an ordinary list of
`True`/`False`.
