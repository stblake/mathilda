# Fourier

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Fourier[list]`**

gives the discrete Fourier transform of a list of complex numbers.

**`Fourier[list, {p1, p2, ...}]`**

returns the specified positions of the discrete Fourier transform.

**`Exp[2 Pi I b (r-1)(s-1)/n], with {a,b} set by the FourierParameters`**

<details>
<summary>Notes</summary>

The transform of a length-n list u is v\[s\] = 1/n^((1-a)/2) Sum\_r u\[r\] option (default {0,1}; {-1,1} data analysis, {1,-1} signal processing). Exact input is first numericalised with N; the list may be a nested rectangular array for a multidimensional transform. Symbolic input yields the exact transform in terms of roots of unity.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Fourier[{1, 1, 2, 2, 1, 1, 0, 0}]
Out[1]= {2.82843, -0.5 + 1.20711*I, 0.0, 0.5 - 0.207107*I, 0.0, 0.5 + 0.207107*I, 0.0, -0.5 - 1.20711*I}

In[2]:= Abs[Fourier[{1, 2, 3, 4, 5, 6}]]^2
Out[2]= {73.5, 6.0, 2.0, 1.5, 2.0, 6.0}

In[3]:= Fourier[{a, b, c, d}]
Out[3]= {1/2 (a + b + c + d), 1/2 (a + I b - c - I d), 1/2 (a - b + c - d), 1/2 (a - I b - c + I d)}
```

### Options (1)

```mathematica
In[4]:= Fourier[{1, 0, 1, 0, 0, 1, 0, 0, 0, 1}, FourierParameters -> {-1, 1}][[1]]
Out[4]= 0.4
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Fourier 1200000 (mixed radix) | 12.7 s | 9.68 s | 8.53 s |
| ListConvolve 100000 x 2048 | 9.71 s | 1.1 s | 10.2 s |
| ListCorrelate 100000 x 2048 | 9.68 s | 1.1 s | 10.4 s |
| Fourier 262143 (awkward size) | 5.91 s | 5.13 s | 4.42 s |
| Fourier 2^18 (262144) | 4.46 s | 2.7 s | 2.35 s |
| InverseFourier 2^18 | 3.37 s | 2.8 s | 2.39 s |

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/), [List](../../other-advanced/List/), [InverseFourier](../../fourier-transforms/InverseFourier/), [FourierDCT](../../fourier-transforms/FourierDCT/), [FourierDST](../../fourier-transforms/FourierDST/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)
