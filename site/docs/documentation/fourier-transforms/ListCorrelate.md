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

## Examples (3)

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

**Attributes:** `Protected`.

## References

**See also:** [ListConvolve](../../fourier-transforms/ListConvolve/), [Times](../../arithmetic/Times/), [Plus](../../arithmetic/Plus/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/fourier-transforms.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/fourier-transforms.md)
- Tests: [`tests/test_compile_linalg.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_linalg.c)
- Tests: [`tests/test_convolutions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_convolutions.c)
- Tests: [`tests/test_correlations.c`](https://github.com/stblake/mathilda/blob/main/tests/test_correlations.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
