# FourierDST

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FourierDST[list]`**

gives the Fourier discrete sine transform (type II) of list.

**`FourierDST[list, m]`**

gives the type-m transform, m one of 1..4 or "I".."IV".

<details>
<summary>Notes</summary>

The four real orthonormal types are self/pair-inverse: I and IV invert themselves; II and III invert each other. Exact input is numericalised with N first; list may be a rectangular nested array (transformed per axis). Machine and arbitrary-precision (MPFR) input are both supported.

</details>

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= FourierDST[{0, 0, 1, 0, 1}]
Out[1]= {0.58541, -0.262866, -0.0854102, -0.425325, 0.894427}

In[2]:= FourierDST[{0, 0, 1, 0, 0}, "IV"]
Out[2]= {0.447214, 0.447214, -0.447214, -0.447214, 0.447214}
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[FourierDCT](../../fourier-transforms/FourierDCT/), [NDArray](../../linear-algebra/NDArray/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)
