---
source: src/ndarray.c
---
**Algorithm.** `DataType` plays two roles, both wired in `ndarray.c`.

As a **function**, `builtin_datatype` reads the element type of a packed NDArray.
With one argument it checks `is_ndarray(arg)`; on a hit it returns
`ndt_to_string(arg->data.ndarray.dtype)` — one of the strings `"float64"`,
`"float32"`, `"complex64"`, `"complex32"`, or `"bool"`. On anything that is not an
NDArray (an ordinary `List`, a scalar) it returns `NULL`, so `DataType[...]` stays
symbolic rather than guessing a type.

As an **option keyword**, `DataType` is the name on the left of the
`DataType -> "float32"` rule that `NDArray[list, DataType -> "..."]` reads when
packing, and it is the single entry in `Options[NDArray]` (default
`DataType -> "float64"`). It is registered `Protected`; its docstring lives
centrally in `info.c`.

**Data structures.** The dtype is a small enum tag stored in the NDArray's
`ndarray` struct (it selects the buffer's element width and the BLAS s/d/c/z
precision). The function emits a freshly allocated `EXPR_STRING`.

**Complexity / limits.** `O(1)` — a tag read, no buffer traversal. Only the packed
`NDArray[...]` surface carries a dtype; a plain list has none, which is why
`DataType` declines on one.
