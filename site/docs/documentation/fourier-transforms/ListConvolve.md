# ListConvolve

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ListConvolve[ker, list]`**

forms the convolution Sum\_r ker\[r\] list\[s-r\] of ker with list.

**`ListConvolve[ker, list, k]`**

aligns element k of ker with each element of list (cyclic).

**`ListConvolve[ker, list, {kL, kR}]`**

sets the overhang: {-1,1} none (default), {1,1}/{-1,-1} maximal at one end, {1,-1} maximal at both.

**`ListConvolve[ker, list, klist, padding]`**

pads list at each end with p, cyclic repetitions of {p1,p2,...}, the list itself (default), or {} for no padding.

**`ListConvolve[ker, list, klist, padding, g, h]`**

uses g in place of Times and h in place of Plus.

**`ListConvolve[ker, list, klist, padding, g, h, lev]`**

works at level lev. ker and list may be multidimensional. Large numeric input uses an FFT; exact and symbolic input is computed directly.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= ListConvolve[{x, y}, {a, b, c, d, e, f}]
Out[1]= {b x + a y, c x + b y, d x + c y, e x + d y, f x + e y}

In[2]:= ListConvolve[{{1, 1}, {1, 1}}, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 1]
Out[2]= {{20, 18, 22}, {14, 12, 16}, {26, 24, 28}}

In[3]:= ListCorrelate[{x, y}, {a, b, c, d, e, f}]
Out[3]= {a x + b y, b x + c y, c x + d y, d x + e y, e x + f y}
```

### Applications (4)

No overhang: output is shorter than the list

```mathematica
In[4]:= ListConvolve[{1, 1}, {1, 2, 3, 4}]
Out[4]= {3, 5, 7}
```

A difference kernel gives successive differences

```mathematica
In[5]:= ListConvolve[{1, -1}, {1, 2, 3, 4}]
Out[5]= {1, 1, 1}
```

Symbolic data is handled by the direct engine

```mathematica
In[6]:= ListConvolve[{x, y}, {a, b, c}]
Out[6]= {b x + a y, c x + b y}
```

Maximal overhang at both ends, cyclic wrap

```mathematica
In[7]:= ListConvolve[{1, 1, 1}, {1, 2, 3, 4, 5}, {1, -1}]
Out[7]= {10, 8, 6, 9, 12, 10, 8}
```

## Algorithm

Mathilda — ListConvolve / ListCorrelate.

See convolutions.h for the high-level description.

With kernel K_r and list a_s:

```text
  ListConvolve  computes  Sum_r K_r a_{s-r}
  ListCorrelate computes  Sum_r K_r a_{s+r}
```

over the alignment window fixed by the overhang parameters {kL, kR}.

Alignment (derived to match the Wolfram Language). Per axis, normalise kL,kR to positive kernel indices KL,KR in 1..m (a negative k maps to m+1+k), with a single integer k meaning {k,k}. For output element t (1-based, 1..L) and kernel element r (1-based, 1..m) the list index touched is

```text
  correlate: j = (t - KL) + r,   L = n + KL - KR
  convolve:  j = (t + KL) - r,   L = n - KL + KR
```

result[t] = h( g(K_r, listval(j)) for r=1..m ), default g=Times, h=Plus, with the h-arguments ordered by ascending j. listval(j) resolves out-of-range j via the padding: cyclic list (default), a constant, a cyclic pad list, or empty (the missing list factor is dropped, giving a single-argument g term).

A fully general direct engine handles every case (symbolic / exact / numeric, every padding and overhang, generalized g/h, and n dimensions). For large numeric inputs with the default Times/Plus a separable FFT fast path is used instead — FFTW for machine precision and the MPFR FFT for arbitrary precision, in both 1-D and n-D (see fourier.h for the shared primitives). The fast path materialises the padded list over the exact index window it needs, which reduces every padding mode to a plain linear convolution computed by a zero-padded FFT product; it then slices out the L outputs.

## Implementation notes

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

**Attributes:** `Protected`.

## References

**See also:** [ListCorrelate](../../fourier-transforms/ListCorrelate/), [Times](../../arithmetic/Times/), [Plus](../../arithmetic/Plus/)

- A. V. Oppenheim and R. W. Schafer, *Discrete-Time Signal Processing*, 3rd ed. (Prentice Hall, 2009), ch. 8 (the circular-convolution / DFT-product theorem).
- M. Frigo and S. G. Johnson, *The Design and Implementation of FFTW3*, Proc. IEEE **93** (2005) 216-231.
- Source: [`src/convolutions.c`](https://github.com/stblake/mathilda/blob/main/src/convolutions.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)
- Tests: [`tests/test_convolutions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_convolutions.c)
- Tests: [`tests/test_correlations.c`](https://github.com/stblake/mathilda/blob/main/tests/test_correlations.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)

## Notes & additional examples

### Notes

`ListConvolve[ker, list]` forms `Sum_r ker[r] list[s-r]` over the alignment window set by the
overhang. With no overhang (the default `{-1, 1}`) the output has length
`Length[list] - Length[ker] + 1`.

The third argument sets the overhang. A single integer `k` aligns element `k` of the kernel
with each list element cyclically; `{kL, kR}` places `ker[[kL]]` against the first output and
`ker[[kR]]` against the last. Common settings: `{-1, 1}` none (default), `{1, 1}`/`{-1, -1}`
maximal at one end, `{1, -1}` maximal at both. A fourth argument gives the padding — a constant,
a cyclic pad list, the list itself (the default, so out-of-range indices wrap), or `{}` for no
padding (the missing list factor is simply dropped). A fifth and sixth replace `Times`/`Plus`
with your own `g`/`h`, and a seventh sets the level.

Data may be symbolic, exact, machine, or arbitrary-precision, and kernel and list may be
multidimensional (the overhang broadcasts over axes, or `{{kL…},{kR…}}` sets each axis). Exact
input stays exact. Large numeric input with the default `Times`/`Plus` is computed by an FFT;
packed kernels and lists are read from their buffers and give a packed result. For
one-dimensional data `ListCorrelate[ker, list]` equals `ListConvolve[Reverse[ker], list]`.
