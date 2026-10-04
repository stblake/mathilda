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

## Examples (9)

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

### Applications (5)

The forward transform of a short real list is complex

```mathematica
In[5]:= Fourier[{1, 2, 3, 4}]
Out[5]= {5.0, -1.0 - 1.0*I, -1.0, -1.0 + 1.0*I}
```

A constant list puts all its energy in the zero-frequency term

```mathematica
In[6]:= Fourier[{1, 1, 1, 1}]
Out[6]= {2.0, 0.0, 0.0, 0.0}
```

A symbolic list gives the exact transform in roots of unity

```mathematica
In[7]:= Fourier[{a, b}]
Out[7]= {(a + b)/Sqrt[2], (a - b)/Sqrt[2]}
```

The data-analysis convention scales by 1/n

```mathematica
In[8]:= Fourier[{1, 2, 3, 4}, FourierParameters -> {-1, 1}]
Out[8]= {2.5, -0.5 - 0.5*I, -0.5, -0.5 + 0.5*I}
```

InverseFourier inverts Fourier

```mathematica
In[9]:= InverseFourier[Fourier[{1, 2, 3, 4}]]
Out[9]= {1.0, 2.0, 3.0, 4.0}
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

**Algorithm.** `builtin_fourier` dispatches through `fourier_core` → `fourier_compute`.
With `FourierParameters -> {a, b}` (default `{0, 1}`) the transform is factored as a
*standard* unnormalised DFT `F[u]_k = Sum_j u_j Exp[+2 Pi I jk/n]` (sign `+1` for
`Fourier`, `-1` for `InverseFourier`), evaluated at the gathered index `(b k) mod n`, then
scaled once by `N^{-(1-a)/2}` (`N` the total element count). The `b`-gather is the identity
when `b = 1` and folds the whole `FourierParameters` `b` — the non-invertible `|b| > 1`
cases included — onto one transform; for a rectangular nested array the transform and the
gather run independently per axis (they are separable). `classify` picks one of three
regimes from the leaves:

1. **machine** (all numeric, no MPFR) — `fourier_fft_machine` runs an in-place FFTW nD plan
   (`fftw_plan_dft`, `FFTW_ESTIMATE`, built per call on the actual buffer; `FFTW_BACKWARD`
   realises the `+2 Pi I` convention) or, without `USE_FFTW`, a naive `O(n^2)` separable
   transform. `machine_build_ndarray` then runs the real-collapse: if
   `max|im| <= 16 N eps * max(|re|,|im|)` the imaginary part is dropped and the result is a
   real (`NDT_FLOAT64`) array.
2. **arbitrary** (any MPFR leaf) — `fourier_fft_mpfr` applies a radix-2 Cooley–Tukey FFT per
   axis for power-of-two lengths and a Bluestein chirp-z transform otherwise, at a working
   precision of `target_bits` plus a guard that grows with `log2 N`; imaginary roundoff
   below `|value| 2^-target_bits` is chopped so a numerically real coefficient prints real.
3. **symbolic** (any non-numeric leaf) — `symbolic_path` builds each output as a `Plus` of
   `u_r Exp[2 Pi I b coeff / n]` roots of unity and hands it to the evaluator, so
   `Fourier[{a, b}]` simplifies to `{(a+b)/Sqrt[2], (a-b)/Sqrt[2]}`.

A visible `NDArray[...]` — always machine numeric — takes
`machine_path_ndarray`, reading the packed buffer directly (`nd_gather_interleaved`,
specialised per dtype; `NDT_COMPLEX64` is one `memcpy`) and returning an `NDArray`; a packed
`List` delists so the surface value is a plain `List`. The two-argument position form
`Fourier[list, {p…}]` computes the full transform and returns `Extract[full, {p…}]`.

**Data structures.** An interleaved `(re, im)` `double` buffer of length `2N` on the machine
path; an `ncpx` (MPFR-complex) buffer on the arbitrary path; `Expr` trees on the symbolic
path. Shapes are carried as row-major `int64` `dims`/`strides` (rank up to
`FOURIER_MAX_RANK = NDARRAY_MAX_RANK`).

**Complexity / limits.** `O(N log N)` on FFTW and on the arbitrary path (Bluestein pads to
the next power of two); `O(N n)` per axis for the naive fallback and the symbolic path. `b`
must be a machine integer for the numeric regimes (a non-integer `b` on numeric data stays
unevaluated); the symbolic regime accepts any exact `b`. `Fourier` compiles over a rank-1
machine array via `fourier_compile`, which always builds an `NDT_COMPLEX64` result
(`force_complex`) because a compiled register carries a static element type and the
interpreter's data-dependent collapse-to-real cannot be expressed there (values identical).

**Attributes:** `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/), [List](../../other-advanced/List/), [InverseFourier](../../fourier-transforms/InverseFourier/), [FourierDCT](../../fourier-transforms/FourierDCT/), [FourierDST](../../fourier-transforms/FourierDST/)

- M. Frigo and S. G. Johnson, *The Design and Implementation of FFTW3*, Proc. IEEE **93** (2005) 216-231.
- J. W. Cooley and J. W. Tukey, *An algorithm for the machine calculation of complex Fourier series*, Math. Comp. **19** (1965) 297-301.
- L. I. Bluestein, *A linear filtering approach to the computation of discrete Fourier transform*, IEEE Trans. Audio Electroacoust. **18** (1970) 451-455.
- Source: [`src/fourier.c`](https://github.com/stblake/mathilda/blob/main/src/fourier.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)
- Tests: [`tests/test_parallelmixedtower.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parallelmixedtower.c)

## Notes & additional examples

### Notes

The zero-frequency component is at **position 1**, not the middle of the list.

The second argument of `FourierParameters -> {a, b}` is the convention. The default `{0, 1}`
is the symmetric `1/Sqrt[n]` normalisation (so `Fourier` and `InverseFourier` are unitary and
mutually inverse); `{-1, 1}` is the data-analysis convention (`1/n` forward, `1` inverse) and
`{1, -1}` the signal-processing one (no normalisation forward). `b` folds in a frequency
re-gather and must be a machine integer for numeric input; a non-integer `b` on numeric data
leaves the call unevaluated.

Input decides the regime automatically. Exact numeric data is numericalised with `N` and runs
a machine-precision FFT (FFTW when built with it); MPFR data runs an arbitrary-precision FFT;
and a list with any non-numeric element gives the exact transform in roots of unity, as
`Fourier[{a, b}]` shows. A machine result that is real to within roundoff is collapsed to a
real list.

A list may be a rectangular nested array for a multidimensional transform (applied
independently per axis), and `Fourier[list, {p1, p2, …}]` returns only the named positions. A
visible `NDArray[...]` is transformed on its packed buffer and comes back an `NDArray`;
`Fourier` also lowers inside `Compile[]`, where it always answers a complex array.
