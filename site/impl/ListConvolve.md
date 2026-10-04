---
references:
  - "A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal Processing*, 3rd ed. (Prentice Hall, 2009), ch. 8 (the circular-convolution / DFT-product theorem)."
  - "M. Frigo and S. G. Johnson, *The Design and Implementation of FFTW3*, Proc. IEEE **93** (2005) 216-231."
source: src/convolutions.c
---
**Algorithm.** `builtin_list_convolve` calls the shared `conv_engine(res, CONV_MODE_CONVOLVE)`.
With kernel `K_r` and list `a_s`, the convolution term index per axis is `j = (t + KL) - r`
(correlation uses `j = (t - KL) + r`), where `KL`/`KR` are the overhang parameters `{kL, kR}`
normalised to positive kernel indices `1..m` (`norm_k`; a single integer `k` means `{k, k}`,
and the default is `{-1, 1}` for convolve / `{1, -1}` for correlate). The output length per
axis is `L = n - KL + KR`. `listval_leaf` resolves an out-of-range `j` through the padding —
cyclic list (default), a constant, a cyclic pad list, or empty (the list factor is dropped) —
and `result[t] = h( g(K_r, a_j) )` with `g = Times`, `h = Plus` by default. Kernel and list may
be multidimensional (`CONV_MAX_RANK = 8`); `{kL,kR}` broadcasts over axes and
`{{kL…},{kR…}}` sets each axis independently.

`classify_leaves` chooses the engine. The fully general `conv_direct` (building an `Expr` per
term and evaluating) handles every case — symbolic, exact, a custom `g`/`h`, any padding or
overhang. For inexact numeric data with the default `Times`/`Plus`, `conv_prefer_fft` compares
the *cost* of the direct engine (`prod(L) prod(m)` fused multiply-adds) against three padded
transforms; above the crossover `conv_fft_machine` (FFTW) or `conv_fft_mpfr` (MPFR radix-2 +
Bluestein) materialises the padded window and computes a linear convolution as a zero-padded
FFT product, otherwise `conv_direct_machine` runs a tight `O(L m)` complex sum. That direct
kernel hoists the kernel multi-index decomposition out of the inner loop, special-cases the
interior (no padding) to a strided — at stride ±1, contiguous and vectorisable — dot product,
and advances the output index as an odometer to avoid a per-axis integer division.

**Data structures.** A `ConvSpec` holds the mode, per-axis `kdims`/`ldims`/`Ldims`,
normalised `KL`/`KR`, row-major strides, borrowed kernel/list/pad leaf arrays *or* packed
`NDArray` operands, and the `g`/`h` heads. The numeric path converts every leaf once into a
`ConvNum` of `(re, im)` `double` arrays (or reads an `NDArray` buffer directly); the FFT path
uses interleaved `(re, im)` `double` buffers (`ncpx` for MPFR) over the linear-convolution
length rounded up to a power of two.

**Complexity / limits.** `O(prod(L) prod(m))` direct, `O(P log P)` for the FFT path (`P` the
padded transform size). A packed kernel, packed list, or both are read straight from their
buffers and a packed real result is handed back (`conv_pack_result`); a packed operand that
reaches the general `Expr` engine (symbolic, custom `g`/`h`, a leaf the fast path cannot
represent) degrades via `ndarray_delist_and_reeval` so the answer is the List path's by
construction. `ListConvolve`/`ListCorrelate` lower inside `Compile[]` (two arrays → array).
