# FourierDST

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FourierDST[list]`**

gives the Fourier discrete sine transform (type II) of list.

**`FourierDST[list, m]`**

gives the type-m transform, m one of 1..4 or "I".."IV".

<details>
<summary>Notes</summary>

The four real orthonormal types are self/pair-inverse: I and IV invert themselves; II and III invert each other. Exact input is numericalised with N first; list may be a rectangular nested array (transformed per axis). Machine and arbitrary-precision (MPFR) input are both supported.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= FourierDST[{0, 0, 1, 0, 1}]
Out[1]= {0.58541, -0.262866, -0.0854102, -0.425325, 0.894427}

In[2]:= FourierDST[{0, 0, 1, 0, 0}, "IV"]
Out[2]= {0.447214, 0.447214, -0.447214, -0.447214, 0.447214}
```

### Applications (3)

The default is the type-II sine transform

```mathematica
In[3]:= FourierDST[{1, 2, 3, 4}]
Out[3]= {3.26641, -1.41421, 1.35299, -1.0}
```

Type I; a type is one of 1..4 or "I".."IV"

```mathematica
In[4]:= FourierDST[{1, 2, 3, 4}, 1]
Out[4]= {4.86624, -2.17625, 1.14876, -0.513743}
```

Type I is its own inverse

```mathematica
In[5]:= FourierDST[FourierDST[{1., 2., 3.}, 1], 1]
Out[5]= {1.0, 2.0, 3.0}
```

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [FourierDCT](../../fourier-transforms/FourierDCT/), [NDArray](../../linear-algebra/NDArray/)

- N. Ahmed, T. Natarajan and K. R. Rao, *Discrete Cosine Transform*, IEEE Trans. Comput. **C-23** (1974) 90-93.
- M. Frigo and S. G. Johnson, *The Design and Implementation of FFTW3*, Proc. IEEE **93** (2005) 216-231.
- Source: [`src/fourier.c`](https://github.com/stblake/mathilda/blob/main/src/fourier.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
- Tests: [`tests/test_compile_transforms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_transforms.c)
- Tests: [`tests/test_fourier.c`](https://github.com/stblake/mathilda/blob/main/tests/test_fourier.c)

## Notes & additional examples

### Notes

`FourierDST[list]` is the real orthonormal discrete sine transform, type II by default;
`FourierDST[list, m]` selects a type (`1`..`4`, or `"I"`..`"IV"`). It is the sine-basis dual
of `FourierDCT` and shares its engine.

Inverse pairing is the same as `FourierDCT`: **I and IV are self-inverse, II and III invert
each other**, so `FourierDST[FourierDST[x, 1], 1]` returns `x`. The transform is orthonormal,
so there is no separate inverse scale.

Real input gives a real result; complex input is transformed component-wise. Exact input is
numericalised with `N`; genuinely symbolic input stays unevaluated. A rectangular nested array
is transformed per axis. Machine precision uses FFTW's `RODFT` real-to-real plans when
available (a direct matrix product otherwise); MPFR input takes a direct arbitrary-precision
matrix path. A visible `NDArray[...]` transforms on its buffer, and `FourierDST` lowers inside
`Compile[]` (real → real).
