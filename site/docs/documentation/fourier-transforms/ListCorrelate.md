# ListCorrelate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ListCorrelate[ker, list]`**

forms the correlation Sum\_r ker\[r\] list\[s+r\] of ker with list.

**`ListCorrelate[ker, list, k]`**

aligns element k of ker with each element of list (cyclic).

**`ListCorrelate[ker, list, {kL, kR}]`**

sets the overhang: {1,-1} none (default), {1,1}/{-1,-1} maximal at one end, {-1,1} maximal at both (negated relative to ListConvolve).

**`ListCorrelate[ker, list, klist, padding]`**

pads list as in ListConvolve.

**`ListCorrelate[ker, list, klist, padding, g, h]`**

uses g in place of Times and h in place of Plus.

**`ListCorrelate[ker, list, klist, padding, g, h, lev]`**

works at level lev. Equivalent to ListConvolve\[Reverse\[ker\], list\].

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

No overhang: sliding pair sums

```mathematica
In[4]:= ListCorrelate[{1, 1}, {1, 2, 3, 4}]
Out[4]= {3, 5, 7}
```

A normalised kernel is a moving average

```mathematica
In[5]:= ListCorrelate[{1, 1, 1}/3, {1, 2, 3, 4, 5}]
Out[5]= {2, 3, 4}
```

A single index argument makes the window cyclic

```mathematica
In[6]:= ListCorrelate[{1, 1}, {1, 2, 3, 4}, 1]
Out[6]= {3, 5, 7, 5}
```

Symbolic data, aligned list[s+r] rather than list[s-r]

```mathematica
In[7]:= ListCorrelate[{x, y}, {a, b, c}]
Out[7]= {a x + b y, b x + c y}
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

**Attributes:** `Protected`.

## References

**See also:** [ListConvolve](../../fourier-transforms/ListConvolve/), [Times](../../arithmetic/Times/), [Plus](../../arithmetic/Plus/)

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

`ListCorrelate[ker, list]` forms `Sum_r ker[r] list[s+r]`. It differs from `ListConvolve` only
in the sign of the kernel offset, so for one-dimensional data `ListCorrelate[ker, list]`
equals `ListConvolve[Reverse[ker], list]`. With no overhang (the default `{1, -1}`) the output
has length `Length[list] - Length[ker] + 1`; a normalised constant kernel therefore gives a
moving average.

The argument forms are identical to `ListConvolve` — a single integer `k` (cyclic alignment of
element `k`), a `{kL, kR}` overhang, a padding argument (constant, cyclic pad list, the list
itself by default, or `{}` for none), and the optional `g`/`h` operators and level. The
overhang settings are *negated* relative to `ListConvolve`: `{1, -1}` is none here and `{-1, 1}`
is maximal at both ends.

Data may be symbolic, exact, machine, or arbitrary-precision, and may be multidimensional.
Exact input stays exact; large numeric input with the default `Times`/`Plus` uses an FFT; a
packed kernel or list is read from its buffer and gives a packed result.
