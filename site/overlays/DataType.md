### Worked examples

```mathematica
In[1]:= DataType[NDArray[{1.0, 2.0, 3.0}]]  (* the default packed element type *)
```

```mathematica
In[1]:= DataType[NDArray[{1, 2, 3}, DataType -> "float32"]]  (* DataType is also the option keyword that set it *)
```

```mathematica
In[1]:= DataType[{1, 2, 3}]  (* a plain list has no dtype, so this stays symbolic *)
```

```mathematica
In[1]:= Options[NDArray]  (* DataType is the one NDArray option, defaulting to float64 *)
```

### Notes

`DataType` is both a reader and an option keyword. As a function, `DataType[a]`
gives the element type of a packed `NDArray` `a` as a string — `"float64"`,
`"float32"`, `"complex64"`, `"complex32"` or `"bool"`; on a non-array it declines
and stays symbolic. As an option, `DataType -> "float32"` is how `NDArray[list,
...]` chooses the packed element type, and it is the sole entry of
`Options[NDArray]`. The five type names map onto BLAS's s/d/c/z precisions.
