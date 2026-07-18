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
