# FourierDST

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FourierDST[list]
    gives the Fourier discrete sine transform (type II) of list.
FourierDST[list, m]
    gives the type-m transform, m one of 1..4 or "I".."IV".

The four real orthonormal types are self/pair-inverse: I and IV invert
themselves; II and III invert each other. Exact input is numericalised
with N first; list may be a rectangular nested array (transformed per
axis). Machine and arbitrary-precision (MPFR) input are both supported.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FourierDST[{0, 0, 1, 0, 1}]
Out[1]= {0.58541, -0.262866, -0.0854102, -0.425325, 0.894427}

In[2]:= FourierDST[{0, 0, 1, 0, 0}, "IV"]
Out[2]= {0.447214, 0.447214, -0.447214, -0.447214, 0.447214}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
