# ListQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ListQ[expr] gives True if expr is a list (head List), False otherwise.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

```mathematica
In[1]:= ListQ[{1, 2, 3}]
Out[1]= True

In[2]:= ListQ[5]
Out[2]= False
```

## Implementation notes

`builtin_listq` (`src/list.c`) returns `True`/`False` according to the `is_listq` helper, i.e. whether the argument is an `EXPR_FUNCTION` whose head is the symbol `List`.

**Attributes:** `Protected`.

## References

**See also:** [VectorQ](../../expression-information/VectorQ/), [MatrixQ](../../expression-information/MatrixQ/), [List](../../other-advanced/List/)

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`ListQ` tests only whether the head is `List`; it does not inspect the elements.
