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

## Examples (3)

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

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [NumericQ](../../expression-information/NumericQ/), [Negative](../../expression-information/Negative/), [NonNegative](../../expression-information/NonNegative/), [NonPositive](../../expression-information/NonPositive/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_positive.c`](https://github.com/stblake/mathilda/blob/main/tests/test_positive.c)
