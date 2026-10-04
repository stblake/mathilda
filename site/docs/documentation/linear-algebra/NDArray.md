# NDArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NDArray[nested_list]`**

Packs a rectangular, machine-precision (Integer/Real) nested list into a dense N-dimensional array (numpy ndarray style). Visibly distinct from List: Head, ListQ, and printing never treat an NDArray as a List. Dimensions gives its shape, ArrayDepth its rank, Length its leading-axis length. Builtins that recognize NDArray (Dot, Plus, Times) use a fast C-level path; results that would need a non-machine-precision entry auto-degrade to an ordinary nested List.

**`NDArray[nested_list, DataType -> "float32"]`**

Packs at the given element type: "float64" (default), "float32", "complex64", "complex32", or "bool" (a list of True/False; "Boolean" is accepted too). DataType\[a\] gives an array's type. A ragged (non-rectangular) list is rejected with an NDArray::ragged warning; an empty or non-machine-precision list stays unevaluated.

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= NDArray[{{1, 2}, {3, 4}}]
Out[1]= NDArray[{{1.0, 2.0}, {3.0, 4.0}}]

In[2]:= Dimensions[NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]
Out[2]= {2, 2, 2}

In[3]:= Depth[NDArray[{{1, 2}, {3, 4}}]]
Out[3]= 3

In[4]:= Dot[NDArray[{{1, 2}, {3, 4}}], NDArray[{{5, 6}, {7, 8}}]]
Out[4]= NDArray[{{19.0, 22.0}, {43.0, 50.0}}]

In[5]:= NDArray[{{1, 2}, {3, 4}}] + NDArray[{{5, 6}, {7, 8}}]
Out[5]= NDArray[{{6.0, 8.0}, {10.0, 12.0}}]

In[6]:= NDArrayQ[NDArray[{1, 2, 3}]]
Out[6]= True

In[7]:= NDArray[{{1, x}, {3, 4}}]
Out[7]= NDArray[{{1, x}, {3, 4}}]
```

### Applications (9)

A dense rank-1 machine-precision array

```mathematica
In[8]:= NDArray[{1., 2., 3.}]
Out[8]= NDArray[{1.0, 2.0, 3.0}]
```

The default float64 dtype widens integers to reals

```mathematica
In[9]:= NDArray[{{1, 2}, {3, 4}}]
Out[9]= NDArray[{{1.0, 2.0}, {3.0, 4.0}}]
```

An int64 buffer keeps the entries exact

```mathematica
In[10]:= NDArray[{1, 2, 3}, DataType -> "int64"]
Out[10]= NDArray[{1, 2, 3}]
```

A visible array's Head is NDArray, never List

```mathematica
In[11]:= Head[NDArray[{1., 2., 3.}]]
Out[11]= NDArray
```

Shape read directly off the dims field

```mathematica
In[12]:= Dimensions[NDArray[{{1, 2}, {3, 4}}]]
Out[12]= {2, 2}
```

A visible NDArray is not a List

```mathematica
In[13]:= ListQ[NDArray[{1., 2.}]]
Out[13]= False
```

Listable heads run element-wise on the buffer

```mathematica
In[14]:= Sin[NDArray[{0., 1., 2.}]]
Out[14]= NDArray[{0.0, 0.841471, 0.909297}]
```

Matrix . vector on the dense buffers

```mathematica
In[15]:= Dot[NDArray[{{1., 2.}, {3., 4.}}], NDArray[{1., 1.}]]
Out[15]= NDArray[{3.0, 7.0}]
```

A reduction reads the buffer and returns a scalar

```mathematica
In[16]:= Total[NDArray[{1., 2., 3., 4.}]]
Out[16]= 10.0
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Transpose then Dot (fused?) | 37.1 s | 63.1 s | 18.6 s |
| Partition window 8, offset 1 | 2.39 s | 32.7 s | 5.24 s |
| Transpose 2000x2000 | 1.58 s | 0.736 s | 2.73 s |
| Take rows 1;;1000 of 2000x2000 | 0.215 s | 0.275 s | 0.213 s |
| column slice m[[All, 1]] | 0.004 s | 0.007 s | 0.001 s |
| ArrayReshape 2x10^6 to 1000x2000 | -- | 0.119 s | 0.216 s |

## Implementation notes

**What the object is, and how packing works.** `NDArray` is a first-class, visibly-distinct dense machine-precision array — an N-dimensional (rank 1..`NDARRAY_MAX_RANK` = 64), rectangular, C-order (row-major) buffer, modelled on numpy's `ndarray`. It is a real `EXPR_NDARRAY` node (`expr.h`'s `NDArrayData`) that stores a flat row-major buffer plus rank and dims directly in the `Expr`, not a nested `List[List[...]]` tree, so a recognising builtin reads the buffer instead of paying the flatten/rebuild cost. One node carries two presentations through its `present_as` field: `NDA_HEAD_NDARRAY` — a visible `NDArray[...]` the user typed, whose `Head` is `NDArray`; and `NDA_HEAD_LIST` — a *packed List*, an ordinary `List` the system chose to store densely, whose `Head` stays `List` and which the user sees as a plain list. `builtin_ndarray` builds the visible form through `ndarray_from_nested_list` and is idempotent on an array argument. Automatic packing (`src/pack.c`): `ndbuild_open` packs a produced buffer only when `pack_enabled()` holds and the element count is at least `pack_min_elements()` (`PACK_MIN_ELEMENTS`, overridable by `MATHILDA_PACK_MIN`; disabled wholesale by `MATHILDA_NO_PACK`). The evaluator's **transparency gate** (`evaluate_step` + `pack.c`'s `AWARE` list) then governs visibility: a head on `AWARE` is handed a packed buffer directly, while a head *not* on it has the packed buffer materialised into one boxed `Expr` per element before it runs — so a packed node never sits inside a plain `List` (the no-nesting invariant that keeps the gate's scan `O(argc)`), and a *visible* `NDArray` is never materialised. The ~80 element-wise-kernel heads opt in automatically via `symtab_set_ndarray_*_kernel`; `nd_present_src2` makes a visible `NDArray` dominate a packed `List` in a binary op's result.

**Data structures: the dense buffer + dtype.** `NDType` (declared in `expr.h`) is `NDT_FLOAT64` (the default, value 0 — `double`), `NDT_FLOAT32` (`float`), `NDT_COMPLEX64` / `NDT_COMPLEX32` (interleaved `(re, im)` pairs), `NDT_INT64` (exact machine integers; compiler-internal but reachable from user syntax via `DataType -> "int64"`), and `NDT_BOOL` (one `uint8_t` per element, materialising as `True`/`False`). The `ndt_get` / `ndt_set` pair is the single choke point that widens/narrows between the stored representation and a machine `(double re, double im)` pair, keeping every generic consumer dtype-agnostic; `ndt_get_i` / `ndt_set_i` are the *exact* int64 accessors used wherever the `double` path (exact only to `2^53`) would lose a value, with the `nd_int64_lossy_hit` tripwire guarding the lossy arm. `ndarray_buffer_element_to_expr` is the one place that decides what head an element materialises as (an exact `Integer` for int64, `expr_new_real` for reals, `Complex[re, im]` for complex), so no unpack path can drift. The default float64 dtype widens integer input to reals; `DataType -> "int64"` preserves them.

**Complexity / limits: what stays packed, what materialises.** A chain of `AWARE` heads keeps the buffer the whole way; the first unaware head materialises it (one `Expr` per element, `O(n)`). Several heads are *deliberately* absent from `AWARE`: `List`, `Rule`, `Association`, `Hold`, `HoldForm` enforce the no-nesting invariant, and `ListQ`, `ArrayQ`, `Append`, `LeafCount`, `Simplify` and similar are correct by omission (marking them aware would be a correctness bug — `LeafCount` would count a `10^6`-element array as one node and skew `Simplify`'s metric). The representation is machine-precision only: any entry that would need a symbol, exact/rational, bigint, or MPFR value degrades the whole value to an ordinary nested `List` rather than forcing a lossy conversion (`ndarray_from_nested_list`).

- `Protected`.
- A ragged (non-rectangular) `list` — unequal sublist shapes, or a mix of
  list and non-list siblings — can never form an array, so `NDArray[list]`
  prints a one-line `NDArray::ragged` warning and stays unevaluated. An empty
  list, a non-machine-precision entry (e.g. a symbol), or a non-list argument
  stays unevaluated silently (the symbolic case may become packable after
  further evaluation).
- `Dot[NDArray[a], NDArray[b]]` contracts the trailing axis of `a` with the
  leading axis of `b` over raw doubles for rank <= 2 operands, giving a new
  `NDArray` (or a bare machine `Real` for a vector.vector contraction). Falls
  back to converting through `Normal` and using the generic tensor path for
  higher-rank operands or a rank mismatch; a genuine shape mismatch (inner
  dimensions disagree) prints `Dot::dotsh` and leaves the call unevaluated.
- `NDArray[a] + NDArray[b]` / `NDArray[a] * NDArray[b]` compute elementwise
  `+`/`*` over raw doubles when both operands are `NDArray` values of
  identical shape. When the operands are all `NDArray` values but of
  disagreeing shape, a one-line `NDArray::shape` warning is printed (naming the
  two shapes) and the sum/product is left unevaluated, mirroring `Dot::dotsh`.
  A mixed `NDArray` + scalar/other operand set instead falls through to the
  generic symbolic `Plus`/`Times` path, treating the `NDArray` as an opaque
  term. numpy-style broadcasting (scalar/array, shape-compatible) is not yet
  implemented.
- Because an `NDArray` is purely numeric, combining one with a **symbolic**
  operand (a bare symbol or any non-numeric expression) can never be carried
  out elementwise. `Plus`/`Times`/`Power` print a one-line `NDArray::sym`
  warning and leave the expression unevaluated: `NDArray[{1., 3.}] + a`,
  `c NDArray[{1., 3.}]`, `NDArray[{1., 3.}]^n`. A numeric scalar operand
  (Integer/Real/Rational/Complex) still broadcasts silently and is unaffected.

**Attributes:** `Protected`.

## References

**See also:** [Real](../../other-advanced/Real/), [DataType](../../other-advanced/DataType/), [SameQ](../../comparisons/SameQ/), [List](../../other-advanced/List/), [MatrixQ](../../expression-information/MatrixQ/), [VectorQ](../../expression-information/VectorQ/), [ListQ](../../expression-information/ListQ/), [Head](../../structural-manipulation/Head/)

- Source: [`src/ndarray.c`](https://github.com/stblake/mathilda/blob/main/src/ndarray.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)
- Tests: [`tests/test_assoc_ops.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_ops.c)
- Tests: [`tests/test_bitwise.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bitwise.c)
- Tests: [`tests/test_characteristicpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_characteristicpolynomial.c)

## Notes & additional examples

### Notes

`NDArray[...]` is a first-class, dense, machine-precision array — numpy's
`ndarray`, stored as a flat row-major buffer with a dtype rather than as a nested
`List` tree. Unlike Mathematica's invisible packed arrays, it is always what it
says it is: its `Head` is `NDArray`, so `ListQ` reports `False`. `Dimensions`,
`Length`, element-wise (`Listable`) heads, `Dot`, and reductions such as `Total`
all read the buffer directly.

The dtype defaults to `float64`, which widens integer input to reals; pass
`DataType -> "int64"` to keep machine integers exact, or `"float32"`,
`"complex64"`, `"complex32"`, or `"bool"`. Any value that cannot be a machine
number (a symbol, an exact rational, a bignum) makes the array degrade to an
ordinary `List` rather than lose precision.

Separately from the visible object, Mathilda packs large ordinary lists into the
same dense representation automatically; those *packed lists* still print and
behave as `List`s. A head that has been checked to read a buffer is handed one
directly, and a head that has not gets the buffer materialised into ordinary
expressions first — so a packed list can never hide inside a plain expression.
