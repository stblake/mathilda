# NDArrayQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NDArrayQ[expr]`**

Gives True if expr is an NDArray object, else False.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

A visible NDArray object

```mathematica
In[1]:= NDArrayQ[NDArray[{1, 2, 3}]]
Out[1]= True
```

A big machine-number list packs transparently

```mathematica
In[2]:= NDArrayQ[Range[10000]]
Out[2]= True
```

A small plain list was never packed

```mathematica
In[3]:= NDArrayQ[{1, 2, 3}]
Out[3]= False
```

A scalar is not an array

```mathematica
In[4]:= NDArrayQ[5]
Out[4]= False
```

## Implementation notes

**Algorithm.** `builtin_ndarrayq` (`src/ndarray.c`) is a one-argument predicate: it returns
`True` when `is_ndarray(arg)` holds and `False` otherwise. `is_ndarray` recognises *both*
packed-array surfaces Mathilda uses — a visible `NDArray[...]` object and a transparently
packed `List` (the dense machine-number buffer that `$AutoArrayPacking` builds behind an
ordinary list head) — so `NDArrayQ` is the one test that sees through the transparency gate.

**Data structures.** None beyond the input `Expr`. The function reads the argument's tag and
packed-buffer metadata through `is_ndarray` and constructs a fresh `True`/`False` symbol.
It is registered `Protected` in `ndarray_init`; a non-unary call returns `NULL`. Its sibling
`PackedArrayQ` tests only the packed-`List` surface (`is_packed_list`), so a visible
`NDArray[...]` makes `PackedArrayQ` `False` while `NDArrayQ` is `True`.

**Complexity / limits.** `O(1)` — it inspects the representation, not the contents, so it
never materialises a packed buffer. Because a packed list is otherwise indistinguishable
from a plain `List` (same head, printed form, elements, ordering, matches), `NDArrayQ` is
the intended way to confirm that a value is actually on the buffer fast path: a *small*
plain list that was never packed gives `False`, while `Range[10000]` under the default
`$AutoArrayPacking` gives `True`.

**Attributes:** `Protected`.

## References

- Source: [`src/ndarray.c`](https://github.com/stblake/mathilda/blob/main/src/ndarray.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_diagonal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_diagonal.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)
- Tests: [`tests/test_graph_algos.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_algos.c)

## Notes & additional examples

### Notes

`NDArrayQ[expr]` gives `True` when `expr` is an NDArray object, and `False` otherwise. It is
an `O(1)` representation test that recognises **both** packed-array surfaces: the visible
`NDArray[...]` head and a transparently packed `List` (the dense machine-number buffer
Mathilda builds behind an ordinary list head).

That dual recognition is what makes it the definitive way to see through the packing gate: a
packed list looks exactly like a plain `List` in every other respect, so only `NDArrayQ`
confirms a value is on the buffer fast path — hence `Range[10000]` is `True` while the small
`{1, 2, 3}` is `False`. Its sibling `PackedArrayQ` tests only the packed-`List` surface, so
it would answer `False` on a visible `NDArray[...]`.
