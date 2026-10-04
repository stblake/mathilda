# NoneTrue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NoneTrue[list, test]`**

Gives True if test\[e\] is True for no element e (True for an empty list). Over an association, tests the values.

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

No element is even, so True

```mathematica
In[4]:= NoneTrue[{1, 3, 5}, EvenQ]
Out[4]= True
```

One even element makes it False

```mathematica
In[5]:= NoneTrue[{1, 2, 3}, EvenQ]
Out[5]= False
```

Nothing exceeds 10

```mathematica
In[6]:= NoneTrue[Range[5], # > 10 &]
Out[6]= True
```

## Implementation notes

**Algorithm.** `builtin_none_true` is `all_any_none_true(res, 2)`, the `mode = 2`
case of the quantifier it shares with `AllTrue` (mode 0) and `AnyTrue` (mode 1).
`NoneTrue` is the logical complement of `AnyTrue`: a bool `NDArray` with an
identity-shaped predicate (`TrueQ`/`Identity`) scans the raw byte buffer and
returns `False` if any byte is set, `True` otherwise. The compiled-predicate path
`pred_quantify`, visible-`NDArray` materialisation, and association-over-values
handling are all shared with the siblings.

The fallback walks the elements building and evaluating `test[e]`;
`NoneTrue` short-circuits to `False` on the first `test[e]` that is `True`. A
result that is neither `True` nor `False` leaves the whole call unevaluated
(`return NULL`), and a `Throw` from the test propagates. With no match the result
is `True`, so the empty list gives `True`.

**Data structures.** The argument array is borrowed; each `test[e]` is built,
evaluated, read for `True`/`False`, and freed. The bool-buffer path reads only the
`unsigned char` NDArray buffer.

**Complexity / limits.** `O(n)` tests with early exit; the bool scan is `O(n)`
bytes. An indeterminate predicate value leaves the call unevaluated, matching
Wolfram.

**Attributes:** `Protected`.

## References

**See also:** [AllTrue](../../functional-programming/AllTrue/), [AnyTrue](../../functional-programming/AnyTrue/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_hypergraph.c`](https://github.com/stblake/mathilda/blob/main/tests/test_hypergraph.c)

## Notes & additional examples

### Notes

`NoneTrue[list, test]` is `True` exactly when `test[e]` is `True` for no element
`e` — the negation of `AnyTrue` — and short-circuits to `False` on the first
match. Over an association the values are tested. The empty list gives `True`.

A `test[e]` that is neither `True` nor `False` leaves the whole call unevaluated.
Over a boolean packed array with `TrueQ` or `Identity` as the test, the check is a
single early-exit scan of the raw buffer.
