# SelectFirst

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SelectFirst[list, pred]`**

Gives the first element e of list for which pred\[e\] is True, or Missing\["NotFound"\]. SelectFirst\[list, pred, default\] uses default. Over an association, tests values and returns the first match.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= SelectFirst[{1, 3, 4, 5, 6}, EvenQ]
Out[1]= 4

In[2]:= FirstCase[{1, 2, 3, 4}, _?EvenQ]
Out[2]= 2

In[3]:= SelectFirst[{1, 3, 5}, EvenQ, None]
Out[3]= None
```

### Applications (4)

The first element passing the predicate

```mathematica
In[4]:= SelectFirst[{1, 3, 5, 6, 7}, EvenQ]
Out[4]= 6
```

A pure function as the predicate

```mathematica
In[5]:= SelectFirst[Range[10], # > 5 &]
Out[5]= 6
```

No match, so the supplied default

```mathematica
In[6]:= SelectFirst[{1, 3, 5}, EvenQ, None]
Out[6]= None
```

No match and no default: Missing["NotFound"]

```mathematica
In[7]:= SelectFirst[{1, 3, 5}, EvenQ]
Out[7]= Missing["NotFound"]
```

## Implementation notes

**Algorithm.** `builtin_select_first` returns the first element of a collection
for which a predicate holds, the one-element analogue of `Select`. A
compiled-predicate fast path (`pred_find_first`) scans a packed real buffer and
returns the first passing element directly, or `PRED_NO_MATCH` to fall to the
default/`Missing` branch — cheap when a match comes early, a full scan otherwise.
A visible `NDArray` is materialised and re-dispatched (`ndstruct_delist_repack`);
an association is handled over its values (`assoc_apply_over_values`).

The fallback loops the elements, building and evaluating `pred[e]`, and returns a
copy of the first element whose test is `True`. A `Throw` from the predicate
propagates. If nothing matches, `SelectFirst[list, pred, default]` returns
`default` and the two-argument form returns `Missing["NotFound"]`.

**Data structures.** The argument array is borrowed; each `pred[e]` is a freshly
built, evaluated, then freed `Expr`, and the match is a copy of the stored
element. No result list is accumulated — the scan stops at the first hit.

**Complexity / limits.** `O(1)` when an early element matches, `O(n)` predicate
evaluations when the match is late or absent.

**Attributes:** `Protected`.

## References

**See also:** [FirstCase](../../functional-programming/FirstCase/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_ndarray_selection.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_selection.c)
- Tests: [`tests/test_pred_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_pred_compile.c)

## Notes & additional examples

### Notes

`SelectFirst[list, pred]` returns the first element `e` with `pred[e]` equal to
`True`, scanning left to right and stopping at the first hit — the single-element
companion to `Select`. With no match it returns `Missing["NotFound"]`, or the
`default` given as a third argument. Over an association the values are tested and
the matching value is returned.
