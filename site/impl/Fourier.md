---
references:
  - "M. Frigo and S. G. Johnson, *The Design and Implementation of FFTW3*, Proc. IEEE **93** (2005) 216-231."
  - "J. W. Cooley and J. W. Tukey, *An algorithm for the machine calculation of complex Fourier series*, Math. Comp. **19** (1965) 297-301."
  - "L. I. Bluestein, *A linear filtering approach to the computation of discrete Fourier transform*, IEEE Trans. Audio Electroacoust. **18** (1970) 451-455."
source: src/fourier.c
---
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
