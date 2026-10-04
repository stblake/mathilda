# ToPackedArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ToPackedArray[list] is ToNDArray[list]: it returns list stored as a dense machine-precision buffer. The result is still a List -- same Head, same printed form, same elements -- but NDArrayQ gives True for it. Provided under Mathematica's name for the same operation; see ToNDArray.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= ToPackedArray[{1., 2., 3.}] === ToNDArray[{1., 2., 3.}]
Out[1]= True
```

### Applications (3)

Mathematica's name for ToNDArray

```mathematica
In[2]:= ToPackedArray[{1., 2., 3., 4.}]
Out[2]= {1.0, 2.0, 3.0, 4.0}
```

The result is a packed List

```mathematica
In[3]:= PackedArrayQ[ToPackedArray[{1., 2., 3.}]]
Out[3]= True
```

Literally the same operation

```mathematica
In[4]:= ToPackedArray[{1., 2., 3.}] === ToNDArray[{1., 2., 3.}]
Out[4]= True
```

## Implementation notes

**Algorithm.** `ToPackedArray` is Mathematica's name (`Developer`​`ToPackedArray`) for
`ToNDArray` and is **the same C builtin** (`builtin_tondarray`) registered under a second name
— not a rule that rewrites to `ToNDArray`. Registering it as a genuine alias rather than a
DownValue means it costs no extra evaluation pass, never shows up in traces, and cannot be
shadowed by a user definition on the target. Every form and option is therefore identical to
`ToNDArray`: an optional trailing `DataType -> "..."` is stripped (`pack_take_dtype`), the one
positional argument is packed by `pack_force_coerce` → `pack_build` with no size threshold and
`coerce = true`, so a mixed machine `Integer`/`Real` list widens to a `float64` buffer while
an all-integer list stays `int64`. `ToPackedArray[list] === ToNDArray[list]` by construction.

**Data structures.** Identical to `ToNDArray`: an `EXPR_NDARRAY` with `present_as =
NDA_HEAD_LIST` (the packed-`List` surface — `Head` stays `List`), backed by a dense row-major
buffer sized `ndt_elem_size(dt) * n`.

**Complexity / limits.** `O(n)` (sniff + flatten) and one allocation; ignores the
automatic-packing threshold; rejects complex, ragged, empty, and non-machine input by
returning the list unchanged. See `ToNDArray` for the dtype-coercion rules.

**Attributes:** `Protected`.

## References

**See also:** [ToNDArray](../../packed-arrays/ToNDArray/)

- Source: [`src/pack.c`](https://github.com/stblake/mathilda/blob/main/src/pack.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`ToPackedArray` is the Wolfram Language's name (`Developer`​`ToPackedArray`) for `ToNDArray`,
provided so code written against that name reads across. It is the **same builtin registered
twice**, not a rule that rewrites to `ToNDArray`, so it costs no extra evaluation pass, never
appears in a trace, and cannot be shadowed by a user definition on the other name. Every form
and option is identical: it stores a rectangular machine-number list as a dense buffer
(invisible except to `NDArrayQ`/`PackedArrayQ`), ignores the automatic-packing threshold,
infers `"int64"`/`"float64"`/`"bool"`, widens a mixed `Integer`/`Real` list to `"float64"`,
takes a `DataType -> "..."` override, and returns the list unchanged when it cannot pack. See
`ToNDArray` for the full dtype rules.
