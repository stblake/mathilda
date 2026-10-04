---
references:
  - "N. Ahmed, T. Natarajan and K. R. Rao, *Discrete Cosine Transform*, IEEE Trans. Comput. **C-23** (1974) 90-93."
  - "M. Frigo and S. G. Johnson, *The Design and Implementation of FFTW3*, Proc. IEEE **93** (2005) 216-231."
source: src/fourier.c
---
**Algorithm.** `builtin_fourier_dct` → `dct_core(res, sine = false)` → `dct_compute`.
`FourierDCT[list]` is the type-II transform; `FourierDCT[list, m]` selects type `m`, parsed
by `parse_dct_type` from the integer `1..4` or the strings `"I".."IV"`. Each type is the real
orthonormal (unitary) matrix `v = M u` whose entries are the cosine formulas in the source
header (e.g. DCT-II is `(1/Sqrt[n]) Cos[Pi (j+1/2) i / n]`); the matrix is real and separable,
so it is applied independently per axis and real input gives real output, while complex input
is transformed on its real and imaginary parts separately. `classify` picks the regime:

- **machine** — with `USE_FFTW`, an `O(n log n)` real-to-real plan (`fftw_plan_r2r` with the
  `REDFT00/10/01/11` kind that `dct_fftw_kind` maps type→kind) followed by the per-axis scale
  `dct_axis_scale` that renormalises FFTW's unnormalised output onto the WMA orthonormal
  convention; without FFTW, a direct `O(n^2)` matrix–vector product per axis
  (`dct_apply_axis_machine`, the matrix carrying the full normalisation).
- **arbitrary** (any MPFR leaf) — `dct_matrix_build_mpfr` builds the `n×n` transform matrix at
  a working precision of `target_bits` plus a `log2 N` guard, and applies it directly per axis
  (`dct_arb_path`).
- **symbolic** — genuinely symbolic input stays unevaluated (declines), matching WMA.

A visible `NDArray[...]` is transformed in place on its packed buffer
(`dct_machine_path_ndarray`) and returns an `NDArray`; a packed `List` delists to a plain
`List`. Axis-length validity is enforced by `dct_axis_ok` (DCT-I needs `n >= 2`, which is the
only type that divides by `n-1`).

**Data structures.** Separate `re`/`im` `double` arrays on the machine path (FFTW runs on each
in place); an `n×n` array of `mpfr_t` plus per-line `mpfr_t` scratch on the arbitrary path;
row-major `int64` `dims`/`strides` up to `FOURIER_MAX_RANK`.

**Complexity / limits.** `O(N log N)` with FFTW, `O(N n)` per axis for the matrix fallback and
the MPFR path. Types pair as inverses — I and IV are self-inverse, II and III invert each
other — so `FourierDCT[FourierDCT[x, 2], 3]` returns `x`. `FourierDCT` compiles over a
machine array: the compiler delegates directly to `builtin_fourier_dct`, which already returns
a deterministic real array for real input (no `force_complex` wrapper is needed).
