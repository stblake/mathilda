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

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= InverseFourier[Fourier[{1, 0, 1, 0, 1, 0}]]
Out[1]= {1.0, 0.0, 1.0, 0.0, 1.0, 0.0}
```

### Applications (4)

A constant spectrum inverts to a single spike

```mathematica
In[2]:= InverseFourier[{1, 1, 1, 1}]
Out[2]= {2.0, 0.0, 0.0, 0.0}
```

Only the zero-frequency term, so a constant signal

```mathematica
In[3]:= InverseFourier[{2, 0, 0, 0}]
Out[3]= {1.0, 1.0, 1.0, 1.0}
```

InverseFourier is the exact inverse of Fourier

```mathematica
In[4]:= InverseFourier[Fourier[{1, 2, 3, 4}]]
Out[4]= {1.0, 2.0, 3.0, 4.0}
```

And the composition the other way is identity too

```mathematica
In[5]:= Fourier[InverseFourier[{1., 2., 3., 4.}]]
Out[5]= {1.0, 2.0, 3.0, 4.0}
```

## Implementation notes

**Algorithm.** `builtin_inverse_fourier` is `fourier_core(res, inverse = 1)`: the same engine
as `Fourier`, run with the opposite exponent sign and the inverse normalisation. With
`FourierParameters -> {a, b}` (default `{0, 1}`) the result is
`u_r = N^{-(1+a)/2} Sum_s v_s Exp[-2 Pi I b (r-1)(s-1)/n]`, so the standard transform uses
`base_sign = -1` (FFTW `FFTW_FORWARD`) and the scale exponent is `-(1+a)/2` instead of
`-(1-a)/2`. Everything else is shared with `Fourier`: `classify` selects the machine
(in-place FFTW nD plan, or the naive `O(n^2)` separable fallback), arbitrary (MPFR radix-2
Cooley–Tukey for power-of-two lengths, Bluestein chirp-z otherwise), or symbolic (exact roots
of unity handed to the evaluator) regime; the `b`-gather folds `FourierParameters` `b`;
`machine_build_ndarray` collapses a numerically real result to a real array. By construction
`InverseFourier[Fourier[list]]` returns `list`.

A visible `NDArray[...]` reads its packed buffer directly and returns an `NDArray`; a packed
`List` delists to a plain `List`; `InverseFourier[list, {p…}]` returns
`Extract[InverseFourier[list], {p…}]`.

**Data structures.** Identical to `Fourier`: an interleaved `(re, im)` `double` buffer of
length `2N` (machine), an `ncpx` MPFR-complex buffer (arbitrary), or `Expr` trees (symbolic),
over row-major `int64` `dims`/`strides` up to `FOURIER_MAX_RANK`.

**Complexity / limits.** `O(N log N)` on the FFT paths, `O(N n)` per axis on the naive and
symbolic paths. `b` must be a machine integer for the numeric regimes. `InverseFourier`
compiles over a rank-1 machine array through `inverse_fourier_compile`, which (like
`fourier_compile`) always emits an `NDT_COMPLEX64` array because the compiled register's
element type is static.

**Attributes:** `Protected`.

## References

**See also:** [Fourier](../../fourier-transforms/Fourier/)

- M. Frigo and S. G. Johnson, *The Design and Implementation of FFTW3*, Proc. IEEE **93** (2005) 216-231.
- J. W. Cooley and J. W. Tukey, *An algorithm for the machine calculation of complex Fourier series*, Math. Comp. **19** (1965) 297-301.
- L. I. Bluestein, *A linear filtering approach to the computation of discrete Fourier transform*, IEEE Trans. Audio Electroacoust. **18** (1970) 451-455.
- Source: [`src/fourier.c`](https://github.com/stblake/mathilda/blob/main/src/fourier.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)

## Notes & additional examples

### Notes

`InverseFourier` is `Fourier` with the opposite exponent sign and the inverse normalisation:
with `FourierParameters -> {a, b}` (default `{0, 1}`) it computes
`u[r] = n^{-(1+a)/2} Sum_s v[s] Exp[-2 Pi I b (r-1)(s-1)/n]`. Under the default convention
both transforms are unitary, so `InverseFourier[Fourier[list]] === list` (and the reverse
composition too).

Everything else matches `Fourier`: the same three regimes (machine FFT, arbitrary-precision
MPFR FFT, exact symbolic), the same `FourierParameters` conventions, the same multidimensional
and `Fourier[list, {p1, …}]` position forms, the same `NDArray` buffer fast path, and the same
`Compile[]` lowering (always a complex array). A machine result real to within roundoff
collapses to a real list.
