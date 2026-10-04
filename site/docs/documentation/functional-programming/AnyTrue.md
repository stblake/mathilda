# AnyTrue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AnyTrue[list, test]`**

Gives True if test\[e\] is True for some element e (False for an empty list). Over an association, tests the values.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= AllTrue[{2, 4, 6}, EvenQ]
Out[1]= True

In[2]:= AnyTrue[{1, 3, 4}, EvenQ]
Out[2]= True

In[3]:= NoneTrue[{1, 3, 5}, EvenQ]
Out[3]= True
```

### Applications (3)

No element passes, so False

```mathematica
In[4]:= AnyTrue[{1, 3, 5}, EvenQ]
Out[4]= False
```

One even element is enough for True

```mathematica
In[5]:= AnyTrue[{1, 2, 3}, EvenQ]
Out[5]= True
```

A pure function as the test

```mathematica
In[6]:= AnyTrue[{1, 3, 5}, # > 4 &]
Out[6]= True
```

## Implementation notes

**Algorithm.** `builtin_any_true` is `all_any_none_true(res, 1)`, sharing one
quantifier with `AllTrue` (mode 0) and `NoneTrue` (mode 2). A bool `NDArray` with
an identity-shaped predicate (`TrueQ`/`Identity`) is answered by scanning the raw
byte buffer for any set byte — `np.any`. Otherwise the compiled-predicate path
`pred_quantify` runs, then a visible `NDArray` is materialised and re-dispatched
(`ndstruct_delist_repack`) and an association is quantified over its values.

The fallback walks the elements building and evaluating `test[e]`; `AnyTrue`
short-circuits to `True` on the first `test[e]` that is `True`. A result that is
neither `True` nor `False` makes the whole call stay unevaluated (`return NULL`),
and a `Throw` from the test propagates. If nothing matches, the result is `False`
— so the empty list gives `False`.

**Data structures.** The argument array is borrowed; each `test[e]` is a built,
evaluated, then freed `Expr`. The bool-buffer path reads only the `unsigned char`
NDArray buffer, closing the sign-predicate → quantifier pipeline without
materialising a `True`/`False` symbol per element.

**Complexity / limits.** `O(n)` tests with early exit (`O(1)` when an early
element matches); the bool scan is `O(n)` bytes, no allocation. An indeterminate
predicate value leaves the call unevaluated.

**Attributes:** `Protected`.

## References

**See also:** [AllTrue](../../functional-programming/AllTrue/), [NoneTrue](../../functional-programming/NoneTrue/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_ndarray.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray.c)

## Notes & additional examples

### Notes

`AnyTrue[list, test]` is `True` exactly when `test[e]` is `True` for at least one
element `e`, and short-circuits on the first match. It is the existential
quantifier dual to `AllTrue`; `NoneTrue` is its negation. Over an association the
values are tested. The empty list gives `False`.

A `test[e]` whose value is neither `True` nor `False` leaves the whole call
unevaluated. Over a boolean packed array with `TrueQ` or `Identity` as the test,
the answer is an early-exit scan of the raw byte buffer.
