# Factor shootout — `final`

Times are the best of 3 in-engine runs after a warm-up (the warm-up alone when it exceeds 5 s). **Bold** is the fastest correct engine. `ratio` is Mathilda / best correct engine. CHECK-FAIL = wrong answer (reconstruction failed, or fewer irreducible factors than another verified engine found); its time is discarded.

## cyclotomic

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| cyclotomic-1 | x^120-1 | 91 µs | 264 µs | 4.2 ms | **82 µs** | 47.4 ms | 1.11 | OK |
| cyclotomic-2 | x^360-1 | **158 µs** | 385 µs | 21.0 ms | 253 µs | 381.5 ms | 1.00 | OK |
| cyclotomic-3 | x^840-1 | **460 µs** | 686 µs | 104.0 ms | 605 µs | 3.71 s | 1.00 | OK |
| cyclotomic-4 | x^1260-1 | **475 µs** | 844 µs | 295.0 ms | 803 µs | 8.58 s | 1.00 | OK |
| cyclotomic-5 | x^2310-1 | **1.9 ms** | 2.0 ms | 3.53 s | 13.2 ms | 38.40 s | 1.00 | OK |
| cyclotomic-6 | x^1155+1 | 1.2 ms | **1.1 ms** | 893.0 ms | 12.7 ms | 4.92 s | 1.06 | OK |

## swinnerton-dyer

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| swinnerton-dyer-1 | deg 16 | **689 µs** | 973 µs | 1.2 ms | 2.5 ms | 41.6 ms | 1.00 | OK |
| swinnerton-dyer-2 | deg 32 | **2.1 ms** | 4.2 ms | 5.2 ms | 16.79 s | 767.2 ms | 1.00 | OK |
| swinnerton-dyer-3 | deg 64 | **9.4 ms** | 18.2 ms | 25.0 ms | TIMEOUT | TIMEOUT | 1.00 | OK |
| swinnerton-dyer-4 | SD16(x)*SD16(x+1) | **1.8 ms** | 3.5 ms | 4.4 ms | 51.26 s | 313.2 ms | 1.00 | OK |

## uni-dense

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| uni-dense-1 | 4x deg 8, 4-bit | **1.2 ms** | 1.8 ms | 4.2 ms | 84.8 ms | 71.4 ms | 1.00 | OK |
| uni-dense-2 | 6x deg 15, 8-bit | **4.2 ms** | 8.2 ms | 47.0 ms | 2.15 s | 2.87 s | 1.00 | OK |
| uni-dense-3 | 8x deg 25, 8-bit | **12.6 ms** | 26.3 ms | 135.0 ms | 22.31 s | 25.34 s | 1.00 | OK |
| uni-dense-4 | 5x deg 60, 10-bit | **21.0 ms** | 40.8 ms | 283.0 ms | 72.11 s | 79.39 s | 1.00 | OK |
| uni-dense-5 | 10x deg 40, 16-bit | **51.8 ms** | 117.0 ms | 621.0 ms | 83.22 s | TIMEOUT | 1.00 | OK |

## uni-bigcoef

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| uni-bigcoef-1 | 2x deg 10, 100-bit | **763 µs** | 1.6 ms | 2.7 ms | 13.0 ms | 102.6 ms | 1.00 | OK |
| uni-bigcoef-2 | 3x deg 15, 300-bit | **4.8 ms** | 12.3 ms | 20.0 ms | 147.0 ms | 847.5 ms | 1.00 | OK |
| uni-bigcoef-3 | 3x deg 30, 600-bit | **30.0 ms** | 57.9 ms | 135.0 ms | 1.70 s | 20.66 s | 1.00 | OK |
| uni-bigcoef-4 | 4x deg 40, 1000-bit | **206.3 ms** | 339.5 ms | 610.0 ms | 13.59 s | TIMEOUT | 1.00 | OK |

## uni-multiplicity

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| uni-multiplicity-1 | d3^10 d5^6 | **145 µs** | 306 µs | 900 µs | 3.4 ms | 62.0 ms | 1.00 | OK |
| uni-multiplicity-2 | d2^20 d4^10 d7^5 | **412 µs** | 1.5 ms | 3.6 ms | 18.9 ms | 141.5 ms | 1.00 | OK |
| uni-multiplicity-3 | d5^15 d9^8 d12^4 d3^25 | **1.4 ms** | 6.7 ms | 17.0 ms | 80.1 ms | 564.5 ms | 1.00 | OK |

## uni-sparse

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| uni-sparse-1 | 2x 4-term deg 300 | **87.6 ms** | 526.1 ms | 1.23 s | TIMEOUT | 101.23 s | 1.00 | OK |
| uni-sparse-2 | 3x 3-term deg 600 | **529.2 ms** | 1.41 s | 2.89 s | TIMEOUT | TIMEOUT | 1.00 | OK |
| uni-sparse-3 | 2x 5-term deg 1500 | **2.27 s** | 6.80 s | 12.16 s | TIMEOUT | TIMEOUT | 1.00 | OK |
| uni-sparse-4 | 3x 4-term deg 2500 | **11.94 s** | 21.96 s | TIMEOUT | TIMEOUT | TIMEOUT | 1.00 | OK |

## mv-dense

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| mv-dense-1 | 2v 3x deg 5 | **985 µs** | 2.1 ms | 501.0 ms | 8.0 ms | 1.32 s | 1.00 | OK |
| mv-dense-2 | 2v 4x deg 10 | **7.7 ms** | 19.6 ms | TIMEOUT | 158.5 ms | 24.45 s | 1.00 | OK |
| mv-dense-3 | 3v 3x deg 4 | **1.9 ms** | 4.9 ms | TIMEOUT | 12.4 ms | 1.26 s | 1.00 | OK |
| mv-dense-4 | 3v 4x deg 6 | **9.4 ms** | 37.0 ms | TIMEOUT | 131.5 ms | 20.58 s | 1.00 | OK |
| mv-dense-5 | 4v 3x deg 3 | **2.3 ms** | 5.5 ms | TIMEOUT | 245.7 ms | 2.06 s | 1.00 | OK |

## mv-sparse

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| mv-sparse-1 | 4v 3x 5-term deg 6 | **883 µs** | 2.1 ms | TIMEOUT | 115.1 ms | 1.11 s | 1.00 | OK |
| mv-sparse-2 | 5v 3x 6-term deg 8 | **9.5 ms** | 19.8 ms | TIMEOUT | 22.60 s | 5.70 s | 1.00 | OK |
| mv-sparse-3 | 6v 3x 6-term deg 10 | **4.2 ms** | 10.5 ms | TIMEOUT | TIMEOUT | 13.40 s | 1.00 | OK |
| mv-sparse-4 | 6v 4x 8-term deg 12 | **97.0 ms** | 304.8 ms | TIMEOUT | TIMEOUT | TIMEOUT | 1.00 | OK |

## mv-nonmonic

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| mv-nonmonic-1 | 3v 2x deg 3, lc deg 2 | **898 µs** | 1.2 ms | 4.14 s | 1.1 ms | 101.7 ms | 1.00 | OK |
| mv-nonmonic-2 | 3v 3x deg 4, lc deg 2 | **2.8 ms** | 5.4 ms | TIMEOUT | 35.8 ms | 671.4 ms | 1.00 | OK |
| mv-nonmonic-3 | 4v 3x deg 4, lc deg 3 | **5.1 ms** | 12.0 ms | TIMEOUT | 76.3 ms | 1.47 s | 1.00 | OK |
| mv-nonmonic-4 | 5v 3x deg 3, lc deg 3 | **4.3 ms** | 9.3 ms | TIMEOUT | 295.5 ms | 1.65 s | 1.00 | OK |

## mv-irreducible

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| mv-irreducible-1 | 2v deg 20 | **1.0 ms** | 1.6 ms | 1.17 s | 4.4 ms | 253.2 ms | 1.00 | OK |
| mv-irreducible-2 | 3v deg 8 | **900 µs** | 1.4 ms | TIMEOUT | 1.9 ms | 161.0 ms | 1.00 | OK |
| mv-irreducible-3 | 4v deg 6 | **793 µs** | 1.7 ms | TIMEOUT | 2.4 ms | 197.3 ms | 1.00 | OK |
| mv-irreducible-4 | 5v deg 4 | **430 µs** | 1.0 ms | TIMEOUT | 1.3 ms | 112.8 ms | 1.00 | OK |

## mv-power

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| mv-power-1 | (x+y+z)^6-1 | **616 µs** | 1.0 ms | 174.0 ms | 951 µs | 56.8 ms | 1.00 | OK |
| mv-power-2 | (x+y+z)^10+1 | **1.1 ms** | 2.1 ms | 28.62 s | 2.6 ms | 158.6 ms | 1.00 | OK |
| mv-power-3 | (x+y+z)^15+1 | **2.7 ms** | 6.0 ms | TIMEOUT | 11.1 ms | 1.77 s | 1.00 | OK |
| mv-power-4 | (x+y+z)^20+1 | **3.1 ms** | 13.5 ms | TIMEOUT | 8.3 ms | 1.09 s | 1.00 | OK |
| mv-power-5 | (x+y+z+w)^8-1 | **1.5 ms** | 5.1 ms | TIMEOUT | 3.4 ms | 463.4 ms | 1.00 | OK |

## mv-determinant

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| mv-determinant-1 | Vandermonde 4 | **245 µs** | 965 µs | 2.9 ms | 434 µs | 40.5 ms | 1.00 | OK |
| mv-determinant-2 | Vandermonde 5 | **557 µs** | 2.1 ms | 939.0 ms | 2.9 ms | 229.1 ms | 1.00 | OK |
| mv-determinant-3 | Vandermonde 6 | **1.7 ms** | 18.0 ms | TIMEOUT | 26.9 ms | 1.94 s | 1.00 | OK |
| mv-determinant-4 | Vandermonde 7 | **9.5 ms** | 224.5 ms | TIMEOUT | 316.7 ms | 11.41 s | 1.00 | OK |
| mv-determinant-5 | Cauchy numerator 3 | **429 µs** | 921 µs | 1.8 ms | 434 µs | 64.7 ms | 1.00 | OK |
| mv-determinant-6 | Cauchy numerator 4 | **1.9 ms** | 11.3 ms | 81.0 ms | 17.6 ms | 1.48 s | 1.00 | OK |

## rational

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| rational-1 | uni deg 40 | 3.1 ms | **2.9 ms** | 3.6 ms | 3.8 ms | 128.3 ms | 1.06 | OK |
| rational-2 | uni deg 150 | **17.4 ms** | 25.1 ms | 61.0 ms | 388.8 ms | 564.9 ms | 1.00 | OK |
| rational-3 | uni deg 240 | **31.6 ms** | 59.3 ms | 85.0 ms | 1.56 s | 1.86 s | 1.00 | OK |
| rational-4 | 2v deg 8 | **671 µs** | 4.6 ms | 4.1 ms | 1.4 ms | 145.0 ms | 1.00 | OK |
| rational-5 | 3v deg 6 | **2.4 ms** | 9.2 ms | 24.0 ms | 3.3 ms | 355.2 ms | 1.00 | OK |
| rational-6 | Cauchy det 3 | **1.2 ms** | 15.5 ms | 23.0 ms | 5.0 ms | 762.1 ms | 1.00 | OK |

## modulus-uni

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| modulus-uni-1 | p=2 deg 1250 | 32.6 ms | 175.2 ms | **15.0 ms** | 50.06 s | TIMEOUT | 2.17 | OK |
| modulus-uni-2 | p=3 deg 625 | **12.1 ms** | 61.7 ms | 116.0 ms | 6.89 s | 24.43 s | 1.00 | OK |
| modulus-uni-3 | p=10007 deg 375 | **18.9 ms** | 59.5 ms | 85.0 ms | TIMEOUT | TIMEOUT | 1.00 | OK |
| modulus-uni-4 | p=2^60-93 deg 187 | **12.3 ms** | 47.8 ms | 78.0 ms | SKIP | 33.71 s | 1.00 | OK |

## modulus-mv

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| modulus-mv-1 | p=2 2v 3x deg 6 | **2.4 ms** | 8.3 ms | 12.02 s | TIMEOUT | SKIP | 1.00 | OK |
| modulus-mv-2 | p=3 3v 3x deg 4 | **5.7 ms** | 15.2 ms | TIMEOUT | CHECK-FAIL | SKIP | 1.00 | OK |
| modulus-mv-3 | p=10007 2v 3x deg 8 | **2.1 ms** | 5.7 ms | 7.04 s | CHECK-FAIL | SKIP | 1.00 | OK |
| modulus-mv-4 | p=2^60-93 3v 3x deg 3 | **1.4 ms** | 3.5 ms | TIMEOUT | SKIP | SKIP | 1.00 | OK |

## gaussian

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| gaussian-1 | x^24-1 | 2.5 ms | 31.7 ms | 1.3 ms | **354 µs** | 1.08 s | 7.01 | OK |
| gaussian-2 | x^120-1 | 14.1 ms | 903.0 ms | 8.2 ms | **1.2 ms** | 25.97 s | 11.52 | OK |
| gaussian-3 | (x^2+k) products | 1.8 ms | 31.3 ms | **1.3 ms** | 8.6 ms | 176.2 ms | 1.36 | OK |
| gaussian-4 | x^2+y^2 products | **1.2 ms** | 8.3 ms | SKIP | 1.3 ms | 166.0 ms | 1.00 | OK |
| gaussian-5 | 3v sums of squares | **1.2 ms** | 4.7 ms | SKIP | 1.2 ms | 1.98 s | 1.00 | OK |

## extension

| case | size | mathilda | reference | pari | maxima | sympy | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| extension-1 | Sqrt[2] uni deg 16 | 1.7 ms | 9.7 ms | **900 µs** | 8.1 ms | 134.3 ms | 1.92 | OK |
| extension-2 | Sqrt[2] uni deg 46 | 7.0 ms | 48.5 ms | **6.5 ms** | 215.5 ms | 1.06 s | 1.07 | OK |
| extension-3 | {Sqrt[2],Sqrt[3]} uni deg 16 | 10.6 ms | 42.5 ms | **1.9 ms** | 27.5 ms | 762.5 ms | 5.58 | OK |
| extension-4 | {Sqrt[2],Sqrt[3]} uni deg 36 | 30.8 ms | 92.9 ms | **6.0 ms** | 407.2 ms | 8.78 s | 5.16 | OK |
| extension-5 | Sqrt[2] 2v | **1.2 ms** | 3.6 ms | SKIP | 2.4 ms | 237.5 ms | 1.00 | OK |
| extension-6 | {Sqrt[2],Sqrt[3]} 3v | 7.4 ms | 40.1 ms | SKIP | **5.0 ms** | TIMEOUT | 1.48 | OK |

## Summary

### Per family (Mathilda / best correct engine)

| family | cases | slower than 1.2× | geo-mean ratio | worst ratio | not OK |
|---|---:|---:|---:|---:|---:|
| cyclotomic | 6 | 0 | 1.03 | 1.11 | 0 |
| swinnerton-dyer | 4 | 0 | 1.00 | 1.00 | 0 |
| uni-dense | 5 | 0 | 1.00 | 1.00 | 0 |
| uni-bigcoef | 4 | 0 | 1.00 | 1.00 | 0 |
| uni-multiplicity | 3 | 0 | 1.00 | 1.00 | 0 |
| uni-sparse | 4 | 0 | 1.00 | 1.00 | 0 |
| mv-dense | 5 | 0 | 1.00 | 1.00 | 0 |
| mv-sparse | 4 | 0 | 1.00 | 1.00 | 0 |
| mv-nonmonic | 4 | 0 | 1.00 | 1.00 | 0 |
| mv-irreducible | 4 | 0 | 1.00 | 1.00 | 0 |
| mv-power | 5 | 0 | 1.00 | 1.00 | 0 |
| mv-determinant | 6 | 0 | 1.00 | 1.00 | 0 |
| rational | 6 | 0 | 1.01 | 1.06 | 0 |
| modulus-uni | 4 | 1 | 1.21 | 2.17 | 0 |
| modulus-mv | 4 | 0 | 1.00 | 1.00 | 0 |
| gaussian | 5 | 3 | 2.56 | 11.52 | 0 |
| extension | 6 | 4 | 2.11 | 5.58 | 0 |

### Cases where Mathilda is slower than the best engine by more than 1.2×

| family | case | size | ratio | Mathilda | best |
|---|---|---|---:|---:|---:|
| modulus-uni | modulus-uni-1 | p=2 deg 1250 | 2.17 | 32.6 ms | 15.0 ms |
| gaussian | gaussian-1 | x^24-1 | 7.01 | 2.5 ms | 354 µs |
| gaussian | gaussian-2 | x^120-1 | 11.52 | 14.1 ms | 1.2 ms |
| gaussian | gaussian-3 | (x^2+k) products | 1.36 | 1.8 ms | 1.3 ms |
| extension | extension-1 | Sqrt[2] uni deg 16 | 1.92 | 1.7 ms | 900 µs |
| extension | extension-3 | {Sqrt[2],Sqrt[3]} uni deg 16 | 5.58 | 10.6 ms | 1.9 ms |
| extension | extension-4 | {Sqrt[2],Sqrt[3]} uni deg 36 | 5.16 | 30.8 ms | 6.0 ms |
| extension | extension-6 | {Sqrt[2],Sqrt[3]} 3v | 1.48 | 7.4 ms | 5.0 ms |

### Mathilda CHECK-FAIL / TIMEOUT / UNSUPPORTED / ERROR

None.

### Verified engines that disagree on the factor count

The lower count is under-factored (reported CHECK-FAIL above).

- modulus-mv-2: mathilda 3, reference 3, maxima 1
- modulus-mv-3: mathilda 3, reference 3, pari 3, maxima 2

## Before and after the review (v0.6060 baseline → v0.6071)

Mathilda / best correct engine, geometric mean per family (the other engines' columns are the same runs):

| family | v0.6060 | v0.6071 | v0.6060 not OK | v0.6071 not OK |
|---|---:|---:|---:|---:|
| cyclotomic | 1420 | 1.03 | 0 | 0 |
| swinnerton-dyer | 1.87 | 1.00 | 0 | 0 |
| uni-dense | 6.30 | 1.00 | 0 | 0 |
| uni-bigcoef | 2.41 | 1.00 | 0 | 0 |
| uni-multiplicity | 57.2 | 1.00 | 0 | 0 |
| uni-sparse | 1.02 | 1.00 | 0 | 0 |
| mv-dense | 3.03 | 1.00 | 0 | 0 |
| mv-sparse | 1.86 | 1.00 | 0 | 0 |
| mv-nonmonic | 2.30 | 1.00 | 0 | 0 |
| mv-irreducible | 2.36 | 1.00 | 0 | 0 |
| mv-power | 1.59 | 1.00 | 0 | 0 |
| mv-determinant | 2.22 | 1.00 | 0 | 0 |
| rational | 3.14 | 1.01 | 2 (under-factored) | 0 |
| modulus-uni | 2.13 | 1.21 | 0 | 0 |
| modulus-mv | — | 1.00 | 4 (unsupported) | 0 |
| gaussian | — | 2.56 | 5 (option ignored) | 0 |
| extension | 128 | 2.11 | 2 (univariate only) | 0 |

Remaining gaps, all on correct answers:

- **`modulus-uni-1`** (GF(2), degree 1250): 32.6 ms against PARI's 15.0 ms. The time is FLINT's `nmod_poly_factor`; PARI has dedicated GF(2)[x] arithmetic.
- **Gaussian** (`x^24 - 1`, `x^120 - 1` with `GaussianIntegers -> True`): Maxima's `gfactor` is fastest. Mathilda is 7–12× behind it and 13–64× ahead of the reference.
- **Two-generator extensions** (`{Sqrt[2], Sqrt[3]}`, univariate): PARI's `nffactor` is 5× faster. Mathilda is 3–4× ahead of the reference.

In both of the last two families the cost is the Trager norm on the large irreducible pieces, plus about 150 µs of fixed set-up per piece.
