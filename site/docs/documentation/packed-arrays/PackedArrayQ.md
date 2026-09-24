# PackedArrayQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PackedArrayQ[expr]`**

Gives True if expr is a packed array -- a List stored as a dense machine-precision buffer -- else False. Provided under Mathematica's name; unlike NDArrayQ, a visible NDArray\[...\] object gives False.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= PackedArrayQ[ToPackedArray[{1., 2., 3.}]]
Out[1]= True

In[2]:= PackedArrayQ[NDArray[{1., 2., 3.}]]
Out[2]= False

In[3]:= NDArrayQ[NDArray[{1., 2., 3.}]]
Out[3]= True

In[4]:= PackedArrayQ[{a, b, c}]
Out[4]= False
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [NDArrayQ](../../other-advanced/NDArrayQ/), [AtomQ](../../expression-information/AtomQ/), [ListQ](../../expression-information/ListQ/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/ndarray.c`](https://github.com/stblake/mathilda/blob/main/src/ndarray.c)
- Specification: [`docs/spec/builtins/packed-arrays.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/packed-arrays.md)
