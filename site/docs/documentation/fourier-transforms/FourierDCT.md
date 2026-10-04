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

## Examples (6)

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

### Applications (3)

The default is the type-II cosine transform

```mathematica
In[4]:= FourierDCT[{1, 2, 3, 4}]
Out[4]= {5.0, -1.57716, 0.0, -0.112085}
```

Type I; a type is one of 1..4 or "I".."IV"

```mathematica
In[5]:= FourierDCT[{1, 2, 3, 4}, 1]
Out[5]= {6.12372, -1.63299, 0.0, -0.408248}
```

Type III inverts type II

```mathematica
In[6]:= FourierDCT[FourierDCT[{1., 2., 3., 4.}, 2], 3]
Out[6]= {1.0, 2.0, 3.0, 4.0}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [N](../../arithmetic/N/), [NDArray](../../linear-algebra/NDArray/)

- N. Ahmed, T. Natarajan and K. R. Rao, *Discrete Cosine Transform*, IEEE Trans. Comput. **C-23** (1974) 90-93.
- M. Frigo and S. G. Johnson, *The Design and Implementation of FFTW3*, Proc. IEEE **93** (2005) 216-231.
- Source: [`src/fourier.c`](https://github.com/stblake/mathilda/blob/main/src/fourier.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)

## Notes & additional examples

### Notes

`FourierDCT[list]` is the real orthonormal (unitary) discrete cosine transform, type II by
default. `FourierDCT[list, m]` selects a type: the integer `1`, `2`, `3`, `4` or the string
`"I"`, `"II"`, `"III"`, `"IV"`.

The types pair as inverses: **I and IV are their own inverses**, while **II and III invert
each other**. So `FourierDCT[FourierDCT[x, 2], 3]` and `FourierDCT[FourierDCT[x, 1], 1]` both
return `x` (up to roundoff). Because the transform is orthonormal there is no separate inverse
normalisation to track.

Real input gives a real result; complex input is transformed on its real and imaginary parts.
Exact input is numericalised with `N` first; a genuinely symbolic list stays unevaluated
(unlike `Fourier`, which has a symbolic regime). The transform is separable, so a rectangular
nested array is transformed per axis. Machine precision uses FFTW's real-to-real plans when
available (a direct matrix product otherwise), and MPFR input takes a direct arbitrary-
precision matrix path. A visible `NDArray[...]` is transformed on its buffer; `FourierDCT`
also lowers inside `Compile[]` (real → real).
