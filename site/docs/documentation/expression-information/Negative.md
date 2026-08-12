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

## Examples (3)

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

## Implementation notes

**Attributes:** `Listable`, `Protected`.

## See also

[NumericQ](../../expression-information/NumericQ/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)
- Tests: [`tests/test_negative.c`](https://github.com/stblake/mathilda/blob/main/tests/test_negative.c)
