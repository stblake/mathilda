# Fourier Transforms

Discrete Fourier transforms of lists of complex numbers and rectangular nested
arrays. `Fourier` computes the forward discrete Fourier transform and
`InverseFourier` its inverse. Both are implemented in `src/fourier.{c,h}` and
dispatch automatically across three regimes — an `O(n log n)` machine-precision
transform via [FFTW](https://www.fftw.org/), a hand-rolled arbitrary-precision
(MPFR) FFT, and an exact symbolic transform — sharing a single core that folds
the `FourierParameters` convention onto one standard transform.

These complement the symbolic [`calculus`](calculus.md) and
[`numerical calculus`](numerical-calculus.md) routines: use a Fourier transform
for spectral analysis of sampled data or arrays.

## Fourier

Discrete Fourier transform of a list of complex numbers (or a rectangular
nested array, for a multidimensional transform). `Fourier[list]` returns a list
the same shape as `list`; `Fourier[list, {p1, p2, …}]` returns only the
positions `{p1, …}` (equivalent to `Extract[Fourier[list], {p1, …}]`).
Implemented in `src/fourier.{c,h}`. Attribute: `Protected`.

With the option `FourierParameters -> {a, b}` (default `{0, 1}`), the transform
of a length-`n` list `u` is

```
v[s] = 1/n^((1-a)/2) Sum[u[r] Exp[2 Pi I b (r-1)(s-1)/n], {r, 1, n}]
```

and `N = Length` (the product of the dimensions) sets the normalisation for a
multidimensional array. The zero-frequency term is at position 1. Common
conventions: `{0, 1}` (default, symmetric `1/Sqrt[n]`), `{-1, 1}` (data
analysis, `1/n`), `{1, -1}` (signal processing, no normalisation). `b` must be
a (nonzero) integer for numeric input; `|b|` relatively prime to `n` keeps the
transform invertible.

Three regimes are dispatched automatically:

- **Machine precision** — an `O(n log n)` transform via **FFTW** (linked when
  built with `USE_FFTW`; a naive `O(n^2)` fallback otherwise). Exact numeric
  input is first numericalised with `N`. A result that is numerically real to
  within roundoff collapses to a real list.
- **Arbitrary precision** — a hand-rolled MPFR-complex FFT (radix-2
  Cooley–Tukey for power-of-two lengths, **Bluestein** chirp-z otherwise),
  applied when any element carries MPFR precision (e.g. `N[list, 24]`).
- **Symbolic** — the exact transform built from roots of unity
  `Exp[2 Pi I b (r-1)(s-1)/n]`, returned when the list contains non-numeric
  elements, e.g. `Fourier[{a, b}]` gives `{(a+b)/Sqrt[2], (a-b)/Sqrt[2]}`.

Given an [`NDArray`](linear-algebra.md) — always machine-precision numeric —
the transform reads the packed buffer directly (no `List` round trip) and
returns an `NDArray`, complex or real-collapsed. `Head[Fourier[NDArray[…]]]` is
`NDArray`.

```
In[1]:= Fourier[{1, 1, 2, 2, 1, 1, 0, 0}]
Out[1]= {2.82843, -0.5 + 1.20711 I, 0., 0.5 - 0.207107 I, 0., 0.5 + 0.207107 I, 0., -0.5 - 1.20711 I}

In[2]:= Abs[Fourier[{1, 2, 3, 4, 5, 6}]]^2
Out[2]= {73.5, 6., 2., 1.5, 2., 6.}

In[3]:= Fourier[{a, b, c, d}]
Out[3]= {1/2 (a + b + c + d), 1/2 (a + I b - c - I d), 1/2 (a - b + c - d), 1/2 (a - I b - c + I d)}

In[4]:= Fourier[{1, 0, 1, 0, 0, 1, 0, 0, 0, 1}, FourierParameters -> {-1, 1}][[1]]
Out[4]= 0.4
```

## InverseFourier

Inverse discrete Fourier transform, the inverse of `Fourier`. Same options,
regimes, and multidimensional/position-form support. Implemented in
`src/fourier.{c,h}`. Attribute: `Protected`.

With `FourierParameters -> {a, b}` (default `{0, 1}`),

```
u[r] = 1/n^((1+a)/2) Sum[v[s] Exp[-2 Pi I b (r-1)(s-1)/n], {s, 1, n}]
```

```
In[1]:= InverseFourier[Fourier[{1, 0, 1, 0, 1, 0}]]
Out[1]= {1., 0., 1., 0., 1., 0.}
```

## FourierDCT

Fourier discrete cosine transform of a real (or complex) list, or a rectangular
nested array for a multidimensional transform. `FourierDCT[list]` gives the
type-II transform; `FourierDCT[list, m]` gives type `m`, where `m` is one of the
integers `1..4` or the strings `"I"`, `"II"`, `"III"`, `"IV"`. Implemented in
`src/fourier.{c,h}`. Attribute: `Protected`.

The four real orthonormal (unitary) types, for a length-`n` list `u` giving `v`:

```
DCT-I    Sqrt[2/(n-1)] (u[1]/2 + Sum[u[r] Cos[Pi/(n-1) (r-1)(s-1)], {r,2,n-1}] + (-1)^(s-1) u[n]/2)
DCT-II   1/Sqrt[n] Sum[u[r] Cos[Pi/n (r-1/2)(s-1)], {r,1,n}]
DCT-III  1/Sqrt[n] (u[1] + 2 Sum[u[r] Cos[Pi/n (r-1)(s-1/2)], {r,2,n}])
DCT-IV   Sqrt[2/n] Sum[u[r] Cos[Pi/n (r-1/2)(s-1/2)], {r,1,n}]
```

The inverse transforms for types 1, 2, 3, 4 are types 1, 3, 2, 4 respectively:
DCT-I and DCT-IV are their own inverses, while DCT-II and DCT-III invert each
other. `FourierDCT[list]` is equivalent to `FourierDCT[list, 2]`. Exact input is
first numericalised with `N`; genuinely symbolic input stays unevaluated.

Two numeric regimes are dispatched automatically: an `O(n log n)`
machine-precision transform via **FFTW**'s real-to-real `REDFT` plans (a direct
`O(n^2)` matrix fallback when built without `USE_FFTW`), and a direct `O(n^2)`
MPFR matrix path when any element carries MPFR precision (e.g. `N[list, 24]`).
An [`NDArray`](linear-algebra.md) argument is transformed in place on its packed
buffer and returns an `NDArray`.

```
In[1]:= FourierDCT[{0, 0, 1, 0, 1}]
Out[1]= {0.894427, -0.425325, -0.0854102, -0.262866, 0.58541}

In[2]:= FourierDCT[{1, 0, 0, 1, 2}, 1]
Out[2]= {1.76777, -0.853553, 1.06066, 0.146447, 0.353553}

In[3]:= FourierDCT[{1, 2 I, 3, 4 I}]
Out[3]= {2. + 3. I, -0.112085 - 1.46508 I, -0.707107 + 0.707107 I, 1.57716 - 1.68925 I}
```

## FourierDST

Fourier discrete sine transform, dual to `FourierDCT`. `FourierDST[list]` gives
the type-II transform; `FourierDST[list, m]` gives type `m` (integer `1..4` or
string `"I".."IV"`). Implemented in `src/fourier.{c,h}`. Attribute: `Protected`.

```
DST-I    Sqrt[2/(n+1)] Sum[u[r] Sin[Pi/(n+1) r s], {r,1,n}]
DST-II   1/Sqrt[n] Sum[u[r] Sin[Pi/n (r-1/2) s], {r,1,n}]
DST-III  1/Sqrt[n] (2 Sum[u[r] Sin[Pi/n r (s-1/2)], {r,1,n-1}] + (-1)^(s-1) u[n])
DST-IV   Sqrt[2/n] Sum[u[r] Sin[Pi/n (r-1/2)(s-1/2)], {r,1,n}]
```

Inverse types match `FourierDCT` (I and IV self-inverse; II and III inverse of
each other). The regimes are identical: FFTW's `RODFT` real-to-real plans for
machine precision, a direct MPFR matrix path for arbitrary precision, and an
`NDArray` fast path.

```
In[1]:= FourierDST[{0, 0, 1, 0, 1}]
Out[1]= {0.58541, -0.262866, -0.0854102, -0.425325, 0.894427}

In[2]:= FourierDST[{0, 0, 1, 0, 0}, "IV"]
Out[2]= {0.447214, 0.447214, -0.447214, -0.447214, 0.447214}
```

## ListConvolve / ListCorrelate

Discrete convolution and correlation, implemented in `src/convolutions.{c,h}`
(shared engine) and `src/correlations.{c,h}` (thin `ListCorrelate` front end).
Attribute on both: `Protected`.

`ListConvolve[ker, list]` forms `Sum_r ker[r] list[s-r]`; `ListCorrelate[ker,
list]` forms `Sum_r ker[r] list[s+r]`, over the alignment window fixed by the
overhang parameters. For one-dimensional data `ListCorrelate[ker, list]` equals
`ListConvolve[Reverse[ker], list]`.

Argument forms (identical for both):

| Form | Meaning |
|------|---------|
| `ListConvolve[ker, list]` | no overhang; output length `Length[list]-Length[ker]+1` |
| `…, k` | cyclic; align element `k` of `ker` with each element (≡ `{k,k}`) |
| `…, {kL, kR}` | overhang at each end: first output has `list[[1]] ker[[kL]]`, last has `list[[-1]] ker[[kR]]` |
| `…, klist, p` | pad ends with `p` |
| `…, klist, {p1,p2,…}` | pad with cyclic repetitions of the `pi` |
| `…, klist, list` | treat `list` as cyclic (default) |
| `…, klist, {}` | no padding (the missing list factor is dropped from each term) |
| `…, klist, padding, g, h` | use `g` in place of `Times`, `h` in place of `Plus` |
| `…, klist, padding, g, h, lev` | operate at level `lev` |

Common `{kL, kR}`: `{-1,1}` no overhang (ListConvolve default), `{1,-1}` no
overhang (ListCorrelate default), `{1,1}`/`{-1,-1}` maximal at one end, and the
"both ends" setting (`{1,-1}` convolve, `{-1,1}` correlate). Settings are negated
between the two functions. Kernels and data may be multidimensional; `{kL,kR}`
broadcasts over axes and `{{kL1,…},{kR1,…}}` sets each axis independently.

Data may be symbolic, exact, machine-precision, or arbitrary-precision. A general
direct engine handles every case (all padding/overhang forms, generalized
`g`/`h`, `n` dimensions). For large numeric input with the default `Times`/`Plus`
a separable FFT fast path is used instead — FFTW for machine precision and the
MPFR FFT (radix-2 + Bluestein) for arbitrary precision, in both 1-D and n-D —
which materialises the padded window and computes a linear convolution via a
zero-padded FFT product. Exact (integer/rational) input stays exact through the
direct path.

```
In[1]:= ListConvolve[{x, y}, {a, b, c, d, e, f}]
Out[1]= {b x + a y, c x + b y, d x + c y, e x + d y, f x + e y}

In[2]:= ListConvolve[{x, y, z}, {1, 2, 3, 4, 5, 6}, {1, -1}]
Out[2]= {x + 6 y + 5 z, 2 x + y + 6 z, 3 x + 2 y + z, 4 x + 3 y + 2 z,
         5 x + 4 y + 3 z, 6 x + 5 y + 4 z, x + 6 y + 5 z, 2 x + y + 6 z}

In[3]:= ListConvolve[{{1, 1}, {1, 1}}, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, 1]
Out[3]= {{20, 18, 22}, {14, 12, 16}, {26, 24, 28}}

In[4]:= ListConvolve[{x, y, z}, {1, 2, 3}, 1, {}, f, g]
Out[4]= {g[f[z], f[y], f[x, 1]], g[f[z], f[y, 1], f[x, 2]],
         g[f[z, 1], f[y, 2], f[x, 3]]}

In[5]:= ListCorrelate[{x, y}, {a, b, c, d, e, f}]
Out[5]= {a x + b y, b x + c y, c x + d y, d x + e y, e x + f y}
```
