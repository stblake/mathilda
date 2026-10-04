### Worked examples

```mathematica
In[1]:= ToPackedArray[{1., 2., 3., 4.}]  (* Mathematica's name for ToNDArray *)
```

```mathematica
In[1]:= PackedArrayQ[ToPackedArray[{1., 2., 3.}]]  (* the result is a packed List *)
```

```mathematica
In[1]:= ToPackedArray[{1., 2., 3.}] === ToNDArray[{1., 2., 3.}]  (* literally the same operation *)
```

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
