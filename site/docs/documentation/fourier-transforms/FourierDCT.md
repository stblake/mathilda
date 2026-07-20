# FourierDCT

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FourierDCT[list]
    gives the Fourier discrete cosine transform (type II) of list.
FourierDCT[list, m]
    gives the type-m transform, m one of 1..4 or "I".."IV".

The four real orthonormal types are self/pair-inverse: I and IV invert
themselves; II and III invert each other. Exact input is numericalised
with N first; list may be a rectangular nested array (transformed per
axis). Machine and arbitrary-precision (MPFR) input are both supported.
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
