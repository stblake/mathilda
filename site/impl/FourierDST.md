---
references:
  - "N. Ahmed, T. Natarajan and K. R. Rao, *Discrete Cosine Transform*, IEEE Trans. Comput. **C-23** (1974) 90-93."
  - "M. Frigo and S. G. Johnson, *The Design and Implementation of FFTW3*, Proc. IEEE **93** (2005) 216-231."
source: src/fourier.c
---
**Algorithm.** `builtin_fourier_dst` is `dct_core(res, sine = true)`: the discrete sine
transform shares `dct_compute` with `FourierDCT`, differing only in the matrix entries and the
FFTW kind. `FourierDST[list]` is the type-II transform; `FourierDST[list, m]` selects type
`m` (`1..4` or `"I".."IV"`, via `parse_dct_type`). Each type is the real orthonormal matrix
`v = M u` built from the sine formulas in the source header (e.g. DST-I is
`Sqrt[2/(n+1)] Sin[Pi (j+1)(i+1)/(n+1)]`), with the DST-III last column special-cased to
`(-1)^i`. The matrix is real and separable, so it is applied independently per axis; real
input gives real output, and complex input is transformed on its real and imaginary parts
separately. Regimes mirror `FourierDCT`:

- **machine** — `fftw_plan_r2r` with the `RODFT00/10/01/11` kind (`dct_fftw_kind` for
  `sine = true`) and the matching `dct_axis_scale` onto the orthonormal convention, or a direct
  `O(n^2)` matrix product per axis without `USE_FFTW`.
- **arbitrary** — `dct_matrix_build_mpfr` builds the `n×n` sine matrix at
  `target_bits + O(log2 N)` working bits and applies it per axis (`dct_arb_path`).
- **symbolic** — declines (stays unevaluated), as WMA does.

A visible `NDArray[...]` transforms in place on its packed buffer and returns an `NDArray`; a
packed `List` delists to a plain `List`.

**Data structures.** The same as `FourierDCT`: separate `re`/`im` `double` arrays (machine),
an `n×n` `mpfr_t` matrix with per-line scratch (arbitrary), row-major `int64`
`dims`/`strides`.

**Complexity / limits.** `O(N log N)` with FFTW, `O(N n)` per axis otherwise. Inverse pairing
matches `FourierDCT`: types I and IV are self-inverse, II and III invert each other, so
`FourierDST[FourierDST[x, 1], 1]` returns `x`. `FourierDST` compiles over a machine array by
direct delegation to `builtin_fourier_dst` (real→real, no complex wrapper).
