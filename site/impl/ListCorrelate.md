---
references:
  - "A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal Processing*, 3rd ed. (Prentice Hall, 2009), ch. 8 (the circular-convolution / DFT-product theorem)."
  - "M. Frigo and S. G. Johnson, *The Design and Implementation of FFTW3*, Proc. IEEE **93** (2005) 216-231."
source: src/convolutions.c
---
**Algorithm.** `builtin_list_correlate` is `conv_engine(res, CONV_MODE_CORRELATE)` — the same
engine as `ListConvolve` with the sign of the kernel offset flipped. With kernel `K_r` and
list `a_s`, the term index per axis is `j = (t - KL) + r` (convolution uses `j = (t + KL) - r`)
and the output length is `L = n + KL - KR`, so for one-dimensional data `ListCorrelate[ker,
list]` equals `ListConvolve[Reverse[ker], list]`. The overhang `{kL, kR}` normalises through
`norm_k` (default `{1, -1}` for correlate; a single integer `k` is `{k, k}`), out-of-range
positions resolve through the same padding modes (cyclic list, constant, cyclic pad list, or
empty), and `result[t] = h( g(K_r, a_j) )` with default `g = Times`, `h = Plus`. Kernel and
list may be multidimensional up to `CONV_MAX_RANK = 8`.

Engine selection is shared: the general `conv_direct` covers symbolic / exact / custom-`g`/`h`
/ any-overhang cases; inexact numeric data with the default heads takes the FFT path
(`conv_fft_machine` via FFTW, or `conv_fft_mpfr` radix-2 + Bluestein) when `conv_prefer_fft`
judges the transform cheaper than the `O(L m)` direct sum, else the vectorised
`conv_direct_machine`. In the FFT path correlation is realised by reversing the kernel on
every axis before the forward transform, so the same zero-padded DFT-product machinery serves
both modes.

**Data structures.** Identical to `ListConvolve`: a `ConvSpec` (mode, per-axis
`kdims`/`ldims`/`Ldims`, normalised `KL`/`KR`, strides, leaf arrays or packed `NDArray`
operands, `g`/`h` heads), a once-converted `ConvNum` of `(re, im)` `double` arrays for the
direct numeric path, and interleaved `(re, im)` `double` (or `ncpx` MPFR) buffers over the
power-of-two padded length for the FFT path.

**Complexity / limits.** `O(prod(L) prod(m))` direct, `O(P log P)` FFT. Packed operands are
read from their buffers and a packed real result is returned; a packed operand that reaches
the general engine degrades through `ndarray_delist_and_reeval`. Exact (integer/rational)
input stays exact on the direct path. `ListCorrelate` lowers inside `Compile[]` (two arrays →
array).
