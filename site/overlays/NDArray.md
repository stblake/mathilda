### Worked examples

```mathematica
In[1]:= NDArray[{1., 2., 3.}]  (* a dense rank-1 machine-precision array *)
```

```mathematica
In[1]:= NDArray[{{1, 2}, {3, 4}}]  (* the default float64 dtype widens integers to reals *)
```

```mathematica
In[1]:= NDArray[{1, 2, 3}, DataType -> "int64"]  (* an int64 buffer keeps the entries exact *)
```

```mathematica
In[1]:= Head[NDArray[{1., 2., 3.}]]  (* a visible array's Head is NDArray, never List *)
```

```mathematica
In[1]:= Dimensions[NDArray[{{1, 2}, {3, 4}}]]  (* shape read directly off the dims field *)
```

```mathematica
In[1]:= ListQ[NDArray[{1., 2.}]]  (* a visible NDArray is not a List *)
```

```mathematica
In[1]:= Sin[NDArray[{0., 1., 2.}]]  (* Listable heads run element-wise on the buffer *)
```

```mathematica
In[1]:= Dot[NDArray[{{1., 2.}, {3., 4.}}], NDArray[{1., 1.}]]  (* matrix . vector on the dense buffers *)
```

```mathematica
In[1]:= Total[NDArray[{1., 2., 3., 4.}]]  (* a reduction reads the buffer and returns a scalar *)
```

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
