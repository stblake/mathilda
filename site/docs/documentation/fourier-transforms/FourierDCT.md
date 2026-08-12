# FourierDCT

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FourierDCT[list]`**

gives the Fourier discrete cosine transform (type II) of list.

**`FourierDCT[list, m]`**

gives the type-m transform, m one of 1..4 or "I".."IV".

<details>
<summary>Notes</summary>

The four real orthonormal types are self/pair-inverse: I and IV invert themselves; II and III invert each other. Exact input is numericalised with N first; list may be a rectangular nested array (transformed per axis). Machine and arbitrary-precision (MPFR) input are both supported.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FourierDCT[{0, 0, 1, 0, 1}]
Out[1]= {0.894427, -0.425325, -0.0854102, -0.262866, 0.58541}

In[2]:= FourierDCT[{1, 0, 0, 1, 2}, 1]
Out[2]= {1.76777, -0.853553, 1.06066, 0.146447, 0.353553}

In[3]:= FourierDCT[{1, 2 I, 3, 4 I}]
Out[3]= {2.0 + 3.0*I, -0.112085 - 1.46508*I, -0.707107 + 0.707107*I, 1.57716 - 1.68925*I}
```

## Implementation notes

**Attributes:** `Protected`.

## See also

[N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)
