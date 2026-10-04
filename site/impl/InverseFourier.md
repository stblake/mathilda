---
references:
  - "M. Frigo and S. G. Johnson, *The Design and Implementation of FFTW3*, Proc. IEEE **93** (2005) 216-231."
  - "J. W. Cooley and J. W. Tukey, *An algorithm for the machine calculation of complex Fourier series*, Math. Comp. **19** (1965) 297-301."
  - "L. I. Bluestein, *A linear filtering approach to the computation of discrete Fourier transform*, IEEE Trans. Audio Electroacoust. **18** (1970) 451-455."
source: src/fourier.c
---
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
