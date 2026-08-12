# NDArray

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NDArray[nested_list]`**

Packs a rectangular, machine-precision (Integer/Real) nested list into a dense N-dimensional array (numpy ndarray style). Visibly distinct from List: Head, ListQ, and printing never treat an NDArray as a List. Dimensions gives its shape, ArrayDepth its rank, Length its leading-axis length. Builtins that recognize NDArray (Dot, Plus, Times) use a fast C-level path; results that would need a non-machine-precision entry auto-degrade to an ordinary nested List.

**`NDArray[nested_list, DataType -> "float32"]`**

Packs at the given element type: "float64" (default), "float32", "complex64", "complex32", or "bool" (a list of True/False; "Boolean" is accepted too). DataType\[a\] gives an array's type. A ragged (non-rectangular) list is rejected with an NDArray::ragged warning; an empty or non-machine-precision list stays unevaluated.

## Examples (7)

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

## See also

[DataType](../../other-advanced/DataType/), [SameQ](../../comparisons/SameQ/), [List](../../other-advanced/List/), [MatrixQ](../../expression-information/MatrixQ/), [VectorQ](../../expression-information/VectorQ/), [ListQ](../../expression-information/ListQ/), [Head](../../structural-manipulation/Head/), [ToNDArray](../../packed-arrays/ToNDArray/)

## References

- Source: [`src/ndarray.c`](https://github.com/stblake/mathilda/blob/main/src/ndarray.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
