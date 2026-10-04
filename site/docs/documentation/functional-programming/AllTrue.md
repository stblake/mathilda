# AllTrue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AllTrue[list, test]`**

Gives True if test\[e\] is True for every element e (True for an empty list). Over an association, tests the values.

## Examples (8)

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

### Applications (5)

Every element passes the test

```mathematica
In[4]:= AllTrue[{2, 4, 6}, EvenQ]
Out[4]= True
```

A pure function as the test

```mathematica
In[5]:= AllTrue[{1, 2, 3, 4}, # > 0 &]
Out[5]= True
```

One odd element is enough for False

```mathematica
In[6]:= AllTrue[{2, 3, 4}, EvenQ]
Out[6]= False
```

Over an association the VALUES are tested

```mathematica
In[7]:= AllTrue[<|a -> 2, b -> 4|>, EvenQ]
Out[7]= True
```

Vacuously True on the empty list

```mathematica
In[8]:= AllTrue[{}, EvenQ]
Out[8]= True
```

## Implementation notes

**Algorithm.** `builtin_all_true` is `all_any_none_true(res, 0)`, the `mode = 0`
case of a shared quantifier (mode 1 is `AnyTrue`, mode 2 is `NoneTrue`). Three
paths, cheapest first. A bool `NDArray` with an identity-shaped predicate (`TrueQ`
or `Identity`, which hand a boolean element back unchanged) is a single
early-exit scan of the raw byte buffer — `np.all` — returning `False` on the first
unset byte. Otherwise the compiled-predicate fast path `pred_quantify` is tried,
then a visible `NDArray` is materialised and re-dispatched (`ndstruct_delist_repack`)
and an association is quantified over its values (`assoc_apply_over_values`).

The fallback walks the elements, building and evaluating `test[e]` for each. It
short-circuits — `AllTrue` returns `False` on the first `test[e]` that is `False`
— and if any result is neither `True` nor `False` the whole call is left
unevaluated (`return NULL`), matching Wolfram rather than guessing. An in-flight
`Throw` from the test propagates out immediately. With no short-circuit, the empty
and all-pass cases give `True`.

**Data structures.** The collection's argument array is borrowed; each `test[e]`
is a freshly built and evaluated `Expr` that is read for `True`/`False` and freed.
The bool-buffer path touches only the `unsigned char` NDArray buffer, so a 10⁶-element
`True`/`False` array never boxes a single symbol.

**Complexity / limits.** `O(n)` element tests with early exit; the bool scan is
`O(n)` over bytes with no allocation. An indeterminate predicate value leaves the
call unevaluated, so partially-symbolic input flows through unchanged.

**Attributes:** `Protected`.

## References

**See also:** [AnyTrue](../../functional-programming/AnyTrue/), [NoneTrue](../../functional-programming/NoneTrue/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_dsolve.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve.c)

## Notes & additional examples

### Notes

`AllTrue[list, test]` is `True` exactly when `test[e]` is `True` for every element
`e`, and short-circuits on the first `False`. It is the universal quantifier to
`AnyTrue`'s existential and `NoneTrue`'s negation. Over an association the values
are tested, not the keys. The empty list is vacuously `True`.

If some `test[e]` returns a value that is neither `True` nor `False`, the whole
call is left unevaluated rather than guessed — so a partially-symbolic list flows
through unchanged. Over a boolean packed array with `TrueQ` or `Identity` as the
test, the answer comes from an early-exit scan of the raw buffer.
