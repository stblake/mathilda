### Worked examples

```mathematica
In[1]:= FromNDArray[ToNDArray[{1., 2., 3.}]]  (* undoes the buffer: an ordinary List of elements *)
```

```mathematica
In[1]:= NDArrayQ[FromNDArray[ToNDArray[{1., 2., 3.}]]]  (* no longer packed *)
```

```mathematica
In[1]:= FromNDArray[NDArray[{1., 2.}]]  (* also unpacks a visible NDArray into a nested List *)
```

### Notes

`FromNDArray[expr]` undoes buffer storage. A packed `List` becomes an ordinary `List` of
separate element nodes (so `NDArrayQ` is then `False`), and an `NDArray[...]` becomes the
nested `List` of its entries. Anything else is returned unchanged. It is the inverse of
`ToNDArray`, and does the same job as `Normal` on these two forms.

The value is identical to the buffer's — only the representation changes — but it allocates one
expression node per element, which is exactly the per-element cost packing exists to avoid, so
this is the point at which a large array stops being cheap to hold.
