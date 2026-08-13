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

## Examples (3)

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

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [NumericQ](../../expression-information/NumericQ/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_nonnegative.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nonnegative.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
