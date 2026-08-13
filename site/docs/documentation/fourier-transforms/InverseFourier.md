# InverseFourier

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`InverseFourier[list]`**

gives the inverse discrete Fourier transform of a list of complex numbers.

**`InverseFourier[list, {p1, p2, ...}]`**

returns the specified positions of the inverse discrete Fourier transform.

<details>
<summary>Notes</summary>

The inverse transform of a length-n list v is u\[r\] = 1/n^((1+a)/2) Sum\_s v\[s\] Exp\[-2 Pi I b (r-1)(s-1)/n\], with {a,b} set by the FourierParameters option (default {0,1}). InverseFourier inverts Fourier.

</details>

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= InverseFourier[Fourier[{1, 0, 1, 0, 1, 0}]]
Out[1]= {1.0, 0.0, 1.0, 0.0, 1.0, 0.0}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Fourier](../../fourier-transforms/Fourier/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)
