# Together / Cancel / Apart shootout — `v6076`

Times are the best of 3 in-engine runs after a warm-up (the warm-up alone when it exceeds 5 s). **Bold** is the fastest correct engine. `ratio` is Mathilda / best correct engine. CHECK-FAIL = wrong or incomplete answer (value differs, not reduced, an Apart term not proper or not over one irreducible); its time is discarded.

## tog-lin

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| tog-lin-1 | n=20 | 100 µs | 3.1 ms | 1.6 ms | 491.2 ms | **71 µs** | 1.41 | OK |
| tog-lin-2 | n=50 | **401 µs** | 33.3 ms | 18.3 ms | 10.97 s | 427 µs | 1.00 | OK |
| tog-lin-3 | n=100 | **1.3 ms** | 222.0 ms | 102.3 ms | TIMEOUT | 2.1 ms | 1.00 | OK |
| tog-lin-4 | n=200 | **5.1 ms** | 1.67 s | 819.8 ms | TIMEOUT | 14.4 ms | 1.00 | OK |

## tog-uni

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| tog-uni-1 | n=8 | 111 µs | 1.3 ms | 607 µs | 136.0 ms | **31 µs** | 3.56 | OK |
| tog-uni-2 | n=25 | 558 µs | 22.8 ms | 22.1 ms | 10.04 s | **328 µs** | 1.70 | OK |
| tog-uni-3 | n=60 | **2.3 ms** | 215.6 ms | 140.6 ms | TIMEOUT | 2.4 ms | 1.00 | OK |
| tog-uni-4 | n=120 | **7.1 ms** | 2.59 s | 1.18 s | TIMEOUT | 13.9 ms | 1.00 | OK |

## tog-multi

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| tog-multi-1 | n=4, 3 vars | **520 µs** | 6.9 ms | 2.0 ms | 907.3 ms | 521 µs | 1.00 | OK |
| tog-multi-2 | n=8, 3 vars | **2.6 ms** | 33.1 ms | 15.1 ms | 33.34 s | 6.2 ms | 1.00 | OK |
| tog-multi-3 | n=12, 3 vars | **9.0 ms** | 139.3 ms | 46.5 ms | TIMEOUT | 30.8 ms | 1.00 | OK |
| tog-multi-4 | n=16, 3 vars | **32.1 ms** | 716.5 ms | 135.1 ms | TIMEOUT | 99.0 ms | 1.00 | OK |

## tog-common

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| tog-common-1 | n=5, shared deg 10 | 116 µs | 2.2 ms | 797 µs | 897.6 ms | **111 µs** | 1.05 | OK |
| tog-common-2 | n=10, shared deg 20 | **378 µs** | 6.8 ms | 3.7 ms | 17.43 s | 467 µs | 1.00 | OK |
| tog-common-3 | n=20, shared deg 30 | **1.1 ms** | 26.8 ms | 18.7 ms | TIMEOUT | 2.3 ms | 1.00 | OK |
| tog-common-4 | n=4, 2 vars | **296 µs** | 3.2 ms | 737 µs | 1.07 s | 301 µs | 1.00 | OK |
| tog-common-5 | n=8, 2 vars | 929 µs | 8.8 ms | 2.4 ms | 23.07 s | **909 µs** | 1.02 | OK |

## tog-telescope

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| tog-telescope-1 | 2n=50 | 17 µs | 146 µs | **13 µs** | 3.6 ms | 85 µs | 1.31 | OK |
| tog-telescope-2 | 2n=200 | 16 µs | 140 µs | **12 µs** | 3.7 ms | 346 µs | 1.33 | OK |
| tog-telescope-3 | 2n=600 | 16 µs | 141 µs | **13 µs** | 3.6 ms | 1.0 ms | 1.23 | OK |
| tog-telescope-4 | 2n=20, 2 vars | 21 µs | 214 µs | **18 µs** | 4.4 ms | 45 µs | 1.17 | OK |
| tog-telescope-5 | 2n=80, 2 vars | 22 µs | 208 µs | **18 µs** | 4.5 ms | 183 µs | 1.22 | OK |

## tog-nested

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| tog-nested-1 | depth 8, 1 var | 46 µs | 396 µs | **21 µs** | 41.1 ms | SKIP | 2.19 | OK |
| tog-nested-2 | depth 16, 1 var | 87 µs | 748 µs | **34 µs** | 239.4 ms | SKIP | 2.56 | OK |
| tog-nested-3 | depth 32, 1 var | 193 µs | 1.6 ms | **75 µs** | TIMEOUT | SKIP | 2.57 | OK |
| tog-nested-4 | depth 8, 2 var | 211 µs | 1.6 ms | **156 µs** | 54.1 ms | SKIP | 1.35 | OK |
| tog-nested-5 | depth 14, 2 var | 471 µs | 4.3 ms | **215 µs** | 353.2 ms | SKIP | 2.19 | OK |

## tog-kernel

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| tog-kernel-1 | n=4 | 314 µs | 3.9 ms | **152 µs** | 36.5 ms | SKIP | 2.07 | OK |
| tog-kernel-2 | n=8 | 869 µs | 5.3 ms | **543 µs** | 72.5 ms | SKIP | 1.60 | OK |
| tog-kernel-3 | n=16 | **3.8 ms** | 26.1 ms | 8.4 ms | 752.1 ms | SKIP | 1.00 | OK |

## tog-gauss

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| tog-gauss-1 | n=5 | 300 µs | 448 µs | **201 µs** | 28.6 ms | SKIP | 1.49 | OK |
| tog-gauss-2 | n=9 | 1.2 ms | 1.0 ms | **479 µs** | 139.9 ms | SKIP | 2.53 | OK |
| tog-gauss-3 | n=17 | **606 µs** | 3.2 ms | CHECK-FAIL | 791.9 ms | SKIP | 1.00 | OK |

## tog-mod

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| tog-mod-1 | n=60, p=7 | 214 µs | 5.9 ms | 863 µs | 13.80 s | **213 µs** | 1.01 | OK |
| tog-mod-2 | n=60, p=1000003 | **376 µs** | 145.2 ms | 29.1 ms | 14.16 s | 709 µs | 1.00 | OK |
| tog-mod-3 | n=60, p=2^61-1 | **593 µs** | 211.4 ms | 28.7 ms | 14.61 s | 752 µs | 1.00 | OK |

## can-uni

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| can-uni-1 | gcd deg 50 | **115 µs** | 1.4 ms | 1.6 ms | 79.0 ms | 304 µs | 1.00 | OK |
| can-uni-2 | gcd deg 200 | **281 µs** | 4.1 ms | 5.5 ms | 279.8 ms | 917 µs | 1.00 | OK |
| can-uni-3 | gcd deg 500 | **597 µs** | 11.2 ms | 13.4 ms | 2.16 s | 2.1 ms | 1.00 | OK |
| can-uni-4 | gcd deg 1000 | **1.1 ms** | 18.5 ms | 26.2 ms | 5.24 s | 4.0 ms | 1.00 | OK |

## can-uni-big

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| can-uni-big-1 | gcd deg 20, 200-bit | **218 µs** | 982 µs | 3.4 ms | 32.8 ms | 274 µs | 1.00 | OK |
| can-uni-big-2 | gcd deg 80, 200-bit | **452 µs** | 3.0 ms | 10.6 ms | 89.8 ms | 641 µs | 1.00 | OK |
| can-uni-big-3 | gcd deg 200, 200-bit | **888 µs** | 4.7 ms | 24.4 ms | 223.5 ms | 1.4 ms | 1.00 | OK |

## can-multi-dense

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| can-multi-dense-1 | 3 vars, gcd deg 4 | **495 µs** | 4.8 ms | 1.6 ms | 171.6 ms | 5.0 ms | 1.00 | OK |
| can-multi-dense-2 | 3 vars, gcd deg 6 | **788 µs** | 8.8 ms | 2.6 ms | 331.3 ms | 8.0 ms | 1.00 | OK |
| can-multi-dense-3 | 3 vars, gcd deg 8 | **1.1 ms** | 14.2 ms | 4.1 ms | 691.8 ms | 10.6 ms | 1.00 | OK |

## can-multi-sparse

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| can-multi-sparse-1 | 6 vars, gcd 10 terms | **846 µs** | 5.2 ms | 4.9 ms | 152.8 ms | 2.7 ms | 1.00 | OK |
| can-multi-sparse-2 | 6 vars, gcd 20 terms | **1.9 ms** | 18.0 ms | 28.2 ms | 572.9 ms | TIMEOUT | 1.00 | OK |
| can-multi-sparse-3 | 6 vars, gcd 40 terms | **5.6 ms** | 67.9 ms | 94.2 ms | 2.54 s | TIMEOUT | 1.00 | OK |

## can-coprime

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| can-coprime-1 | deg 100 | 226 µs | 2.7 ms | 4.7 ms | 115.9 ms | **173 µs** | 1.31 | OK |
| can-coprime-2 | deg 400 | **994 µs** | 10.6 ms | 50.4 ms | 598.6 ms | 1.9 ms | 1.00 | OK |
| can-coprime-3 | deg 1000 | **2.3 ms** | 27.2 ms | 323.9 ms | 3.55 s | 7.6 ms | 1.00 | OK |
| can-coprime-4 | 3 vars deg 5 | **177 µs** | 3.2 ms | 418 µs | 76.7 ms | 5.6 ms | 1.00 | OK |
| can-coprime-5 | 3 vars deg 8 | **580 µs** | 9.8 ms | 1.2 ms | 213.8 ms | 377.0 ms | 1.00 | OK |

## can-small

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| can-small-1 | textbook | 12 µs | 71 µs | 7 µs | 2.2 ms | **1 µs** | 8.80 | OK |
| can-small-2 | textbook | 17 µs | 152 µs | 14 µs | 3.9 ms | **2 µs** | 8.56 | OK |
| can-small-3 | textbook | 16 µs | 129 µs | 20 µs | 2.5 ms | **8 µs** | 2.08 | OK |
| can-small-4 | textbook | 16 µs | 150 µs | 13 µs | 3.1 ms | **2 µs** | 8.92 | OK |
| can-small-5 | textbook | 18 µs | 163 µs | 29 µs | 5.2 ms | **7 µs** | 2.46 | OK |
| can-small-6 | textbook | 15 µs | 128 µs | 11 µs | 3.7 ms | **2 µs** | 6.90 | OK |

## can-gauss

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| can-gauss-1 | gcd deg 5 | 212 µs | 1.1 ms | **169 µs** | 31.9 ms | SKIP | 1.25 | OK |
| can-gauss-2 | gcd deg 20 | 679 µs | 5.1 ms | **641 µs** | 71.8 ms | SKIP | 1.06 | OK |
| can-gauss-3 | gcd deg 60 | 2.9 ms | 34.9 ms | **2.8 ms** | 189.6 ms | SKIP | 1.02 | OK |

## can-ext

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| can-ext-1 | cofactor deg 5 | **743 µs** | 2.0 ms | CHECK-FAIL | 23.1 ms | SKIP | 1.00 | OK |
| can-ext-2 | cofactor deg 20 | **2.0 ms** | 626.4 ms | CHECK-FAIL | 83.1 ms | SKIP | 1.00 | OK |
| can-ext-3 | cofactor deg 60 | **5.1 ms** | 5.36 s | CHECK-FAIL | 438.6 ms | SKIP | 1.00 | OK |

## can-mod

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| can-mod-1 | gcd deg 200, p=7 | 204 µs | 4.3 ms | 12.0 ms | 316.5 ms | **108 µs** | 1.89 | OK |
| can-mod-2 | gcd deg 200, p=1000003 | 214 µs | 4.7 ms | 5.3 ms | 325.4 ms | **104 µs** | 2.05 | OK |
| can-mod-3 | gcd deg 200, p=2^61-1 | 410 µs | 4.9 ms | 7.3 ms | 295.5 ms | **256 µs** | 1.60 | OK |

## ap-lin

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| ap-lin-1 | 10 linear | **354 µs** | 631 µs | 820 µs | ERROR | SKIP | 1.00 | OK |
| ap-lin-2 | 30 linear | **2.2 ms** | 6.0 ms | 13.3 ms | ERROR | SKIP | 1.00 | OK |
| ap-lin-3 | 60 linear | **9.5 ms** | 34.6 ms | 61.9 ms | ERROR | SKIP | 1.00 | OK |
| ap-lin-4 | 100 linear | **29.2 ms** | 109.7 ms | 323.2 ms | ERROR | SKIP | 1.00 | OK |

## ap-mult

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| ap-mult-1 | mult 3,3,3 | **119 µs** | 467 µs | 305 µs | ERROR | SKIP | 1.00 | OK |
| ap-mult-2 | mult 10,8,5 | **408 µs** | 4.2 ms | 862 µs | ERROR | SKIP | 1.00 | OK |
| ap-mult-3 | mult 30,20,10 | **2.8 ms** | 32.6 ms | 8.3 ms | ERROR | SKIP | 1.00 | OK |
| ap-mult-4 | mult 60,40,20 | **28.7 ms** | 15.13 s | 162.4 ms | ERROR | SKIP | 1.00 | OK |

## ap-irr

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| ap-irr-1 | 2 irreducibles, mult 1 | **88 µs** | 246 µs | 131 µs | ERROR | SKIP | 1.00 | OK |
| ap-irr-2 | 4 irreducibles, mult 2 | **343 µs** | 2.1 ms | 915 µs | ERROR | SKIP | 1.00 | OK |
| ap-irr-3 | 6 irreducibles, mult 3 | **1.4 ms** | 17.6 ms | 5.7 ms | ERROR | SKIP | 1.00 | OK |

## ap-improper

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| ap-improper-1 | 5 linear, num deg 15 | **179 µs** | 543 µs | 352 µs | ERROR | SKIP | 1.00 | OK |
| ap-improper-2 | 20 linear, num deg 30 | **1.2 ms** | 3.1 ms | 5.4 ms | ERROR | SKIP | 1.00 | OK |
| ap-improper-3 | 50 linear, num deg 60 | **6.4 ms** | 24.5 ms | 41.1 ms | ERROR | SKIP | 1.00 | OK |

## ap-param

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| ap-param-1 | 2 params | 198 µs | 264 µs | **43 µs** | ERROR | SKIP | 4.60 | OK |
| ap-param-2 | 3 params | 486 µs | 802 µs | **104 µs** | ERROR | SKIP | 4.67 | OK |
| ap-param-3 | 4 params | 977 µs | 1.4 ms | **602 µs** | ERROR | SKIP | 1.62 | OK |
| ap-param-4 | 5 params | **1.8 ms** | 3.3 ms | 4.4 ms | ERROR | SKIP | 1.00 | OK |
| ap-param-5 | 6 params | **3.2 ms** | 4.2 ms | 37.0 ms | ERROR | SKIP | 1.00 | OK |
| ap-param-6 | (x+a)^2 (x^2+b) | 499 µs | 2.7 ms | **157 µs** | ERROR | SKIP | 3.18 | OK |

## ap-param-content

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| ap-param-content-1 | 4 linear, param in num/content | 314 µs | 519 µs | **245 µs** | ERROR | SKIP | 1.28 | OK |
| ap-param-content-2 | 10 linear, param in num/content | 1.1 ms | 1.2 ms | **1.1 ms** | ERROR | SKIP | 1.00 | OK |
| ap-param-content-3 | 20 linear, param in num/content | 3.5 ms | **3.0 ms** | 6.0 ms | ERROR | SKIP | 1.14 | OK |

## ap-gauss

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| ap-gauss-1 | 3 gaussian linear | 509 µs | **300 µs** | CHECK-FAIL | ERROR | SKIP | 1.70 | OK |
| ap-gauss-2 | 6 gaussian linear | 996 µs | **804 µs** | CHECK-FAIL | ERROR | SKIP | 1.24 | OK |
| ap-gauss-3 | 10 gaussian linear | **1.8 ms** | 6.0 ms | CHECK-FAIL | ERROR | SKIP | 1.00 | OK |

## ap-mod

| case | size | mathilda | reference | maxima | sympy | pari | ratio | Mathilda |
|---|---|---:|---:|---:|---:|---:|---:|---|
| ap-mod-1 | 6 linear, p=7 | **107 µs** | 312 µs | 123 µs | SKIP | SKIP | 1.00 | OK |
| ap-mod-2 | 6 linear, p=1000003 | 194 µs | 370 µs | **171 µs** | SKIP | SKIP | 1.13 | OK |
| ap-mod-3 | 6 linear, p=2^61-1 | 409 µs | 674 µs | **195 µs** | SKIP | SKIP | 2.10 | OK |

## Summary

### Per family (Mathilda / best correct engine)

| family | cases | slower than 1.2× | geo-mean ratio | worst ratio | not OK |
|---|---:|---:|---:|---:|---:|
| tog-lin | 4 | 1 | 1.09 | 1.41 | 0 |
| tog-uni | 4 | 2 | 1.57 | 3.56 | 0 |
| tog-multi | 4 | 0 | 1.00 | 1.00 | 0 |
| tog-common | 5 | 0 | 1.01 | 1.05 | 0 |
| tog-telescope | 5 | 4 | 1.25 | 1.33 | 0 |
| tog-nested | 5 | 5 | 2.12 | 2.57 | 0 |
| tog-kernel | 3 | 2 | 1.49 | 2.07 | 0 |
| tog-gauss | 3 | 2 | 1.56 | 2.53 | 0 |
| tog-mod | 3 | 0 | 1.00 | 1.01 | 0 |
| can-uni | 4 | 0 | 1.00 | 1.00 | 0 |
| can-uni-big | 3 | 0 | 1.00 | 1.00 | 0 |
| can-multi-dense | 3 | 0 | 1.00 | 1.00 | 0 |
| can-multi-sparse | 3 | 0 | 1.00 | 1.00 | 0 |
| can-coprime | 5 | 1 | 1.05 | 1.31 | 0 |
| can-small | 6 | 6 | 5.36 | 8.92 | 0 |
| can-gauss | 3 | 1 | 1.11 | 1.25 | 0 |
| can-ext | 3 | 0 | 1.00 | 1.00 | 0 |
| can-mod | 3 | 3 | 1.84 | 2.05 | 0 |
| ap-lin | 4 | 0 | 1.00 | 1.00 | 0 |
| ap-mult | 4 | 0 | 1.00 | 1.00 | 0 |
| ap-irr | 3 | 0 | 1.00 | 1.00 | 0 |
| ap-improper | 3 | 0 | 1.00 | 1.00 | 0 |
| ap-param | 6 | 4 | 2.19 | 4.67 | 0 |
| ap-param-content | 3 | 1 | 1.14 | 1.28 | 0 |
| ap-gauss | 3 | 2 | 1.28 | 1.70 | 0 |
| ap-mod | 3 | 1 | 1.34 | 2.10 | 0 |

### Cases where Mathilda is slower than the best engine by more than 1.2×

| family | case | size | ratio | Mathilda | best |
|---|---|---|---:|---:|---:|
| tog-lin | tog-lin-1 | n=20 | 1.41 | 100 µs | 71 µs |
| tog-uni | tog-uni-1 | n=8 | 3.56 | 111 µs | 31 µs |
| tog-uni | tog-uni-2 | n=25 | 1.70 | 558 µs | 328 µs |
| tog-telescope | tog-telescope-1 | 2n=50 | 1.31 | 17 µs | 13 µs |
| tog-telescope | tog-telescope-2 | 2n=200 | 1.33 | 16 µs | 12 µs |
| tog-telescope | tog-telescope-3 | 2n=600 | 1.23 | 16 µs | 13 µs |
| tog-telescope | tog-telescope-5 | 2n=80, 2 vars | 1.22 | 22 µs | 18 µs |
| tog-nested | tog-nested-1 | depth 8, 1 var | 2.19 | 46 µs | 21 µs |
| tog-nested | tog-nested-2 | depth 16, 1 var | 2.56 | 87 µs | 34 µs |
| tog-nested | tog-nested-3 | depth 32, 1 var | 2.57 | 193 µs | 75 µs |
| tog-nested | tog-nested-4 | depth 8, 2 var | 1.35 | 211 µs | 156 µs |
| tog-nested | tog-nested-5 | depth 14, 2 var | 2.19 | 471 µs | 215 µs |
| tog-kernel | tog-kernel-1 | n=4 | 2.07 | 314 µs | 152 µs |
| tog-kernel | tog-kernel-2 | n=8 | 1.60 | 869 µs | 543 µs |
| tog-gauss | tog-gauss-1 | n=5 | 1.49 | 300 µs | 201 µs |
| tog-gauss | tog-gauss-2 | n=9 | 2.53 | 1.2 ms | 479 µs |
| can-coprime | can-coprime-1 | deg 100 | 1.31 | 226 µs | 173 µs |
| can-small | can-small-1 | textbook | 8.80 | 12 µs | 1 µs |
| can-small | can-small-2 | textbook | 8.56 | 17 µs | 2 µs |
| can-small | can-small-3 | textbook | 2.08 | 16 µs | 8 µs |
| can-small | can-small-4 | textbook | 8.92 | 16 µs | 2 µs |
| can-small | can-small-5 | textbook | 2.46 | 18 µs | 7 µs |
| can-small | can-small-6 | textbook | 6.90 | 15 µs | 2 µs |
| can-gauss | can-gauss-1 | gcd deg 5 | 1.25 | 212 µs | 169 µs |
| can-mod | can-mod-1 | gcd deg 200, p=7 | 1.89 | 204 µs | 108 µs |
| can-mod | can-mod-2 | gcd deg 200, p=1000003 | 2.05 | 214 µs | 104 µs |
| can-mod | can-mod-3 | gcd deg 200, p=2^61-1 | 1.60 | 410 µs | 256 µs |
| ap-param | ap-param-1 | 2 params | 4.60 | 198 µs | 43 µs |
| ap-param | ap-param-2 | 3 params | 4.67 | 486 µs | 104 µs |
| ap-param | ap-param-3 | 4 params | 1.62 | 977 µs | 602 µs |
| ap-param | ap-param-6 | (x+a)^2 (x^2+b) | 3.18 | 499 µs | 157 µs |
| ap-param-content | ap-param-content-1 | 4 linear, param in num/content | 1.28 | 314 µs | 245 µs |
| ap-gauss | ap-gauss-1 | 3 gaussian linear | 1.70 | 509 µs | 300 µs |
| ap-gauss | ap-gauss-2 | 6 gaussian linear | 1.24 | 996 µs | 804 µs |
| ap-mod | ap-mod-3 | 6 linear, p=2^61-1 | 2.10 | 409 µs | 195 µs |

### Mathilda CHECK-FAIL / TIMEOUT / UNSUPPORTED / ERROR

None.

## Before and after (v0.6071 → v0.6072 → v0.6076)

The largest case of each family. The other engines' times are the same runs in every column. v0.6076 renders the factored form (v0.6073); v0.6072 rendered the expanded one.

| case | size | v0.6071 | v0.6072 | v0.6076 | best other engine |
|---|---|---:|---:|---:|---|
| tog-lin-4 | n=200 | 11.9 ms | 10.0 ms | 5.1 ms | 14.4 ms (pari) |
| tog-uni-4 | n=120 | 17.7 ms | 12.5 ms | 7.1 ms | 13.9 ms (pari) |
| tog-multi-4 | n=16, 3 vars | 136.3 ms | 91.2 ms | 32.1 ms | 99.0 ms (pari) |
| tog-common-5 | n=8, 2 vars | 1.5 ms | 1.2 ms | 929 µs | 909 µs (pari) |
| tog-telescope-5 | 2n=80, 2 vars | 39 µs | 30 µs | 22 µs | 18 µs (maxima) |
| tog-nested-5 | depth 14, 2 var | 724 µs | 536 µs | 471 µs | 215 µs (maxima) |
| tog-kernel-3 | n=16 | 5.2 ms | 4.3 ms | 3.8 ms | 8.4 ms (maxima) |
| tog-gauss-3 | n=17 | 44.4 ms | 479 µs | 606 µs | 3.2 ms (reference) |
| tog-mod-3 | n=60, p=2^61-1 | UNSUPPORTED | 673 µs | 593 µs | 752 µs (pari) |
| can-uni-4 | gcd deg 1000 | 88.5 ms | 1.2 ms | 1.1 ms | 4.0 ms (pari) |
| can-uni-big-3 | gcd deg 200, 200-bit | 6.7 ms | 936 µs | 888 µs | 1.4 ms (pari) |
| can-multi-dense-3 | 3 vars, gcd deg 8 | 4.6 ms | 1.1 ms | 1.1 ms | 4.1 ms (maxima) |
| can-multi-sparse-3 | 6 vars, gcd 40 terms | 20.6 ms | 5.8 ms | 5.6 ms | 67.9 ms (reference) |
| can-coprime-5 | 3 vars deg 8 | 3.0 ms | 1.0 ms | 580 µs | 1.2 ms (maxima) |
| can-small-6 | textbook | 30 µs | 17 µs | 15 µs | 2 µs (pari) |
| can-gauss-3 | gcd deg 60 | 24.0 ms | 2.9 ms | 2.9 ms | 2.8 ms (maxima) |
| can-ext-3 | cofactor deg 60 | 3.79 s | 5.2 ms | 5.1 ms | 438.6 ms (sympy) |
| can-mod-3 | gcd deg 200, p=2^61-1 | UNSUPPORTED | 439 µs | 410 µs | 256 µs (pari) |
| ap-lin-4 | 100 linear | 40.2 ms | 31.0 ms | 29.2 ms | 109.7 ms (reference) |
| ap-mult-4 | mult 60,40,20 | 39.5 ms | 32.0 ms | 28.7 ms | 162.4 ms (maxima) |
| ap-irr-3 | 6 irreducibles, mult 3 | 2.2 ms | 1.4 ms | 1.4 ms | 5.7 ms (maxima) |
| ap-improper-3 | 50 linear, num deg 60 | 41.3 ms | 7.4 ms | 6.4 ms | 24.5 ms (reference) |
| ap-param-6 | (x+a)^2 (x^2+b) | 4.1 ms | 487 µs | 499 µs | 157 µs (maxima) |
| ap-param-content-3 | 20 linear, param in num/content | 64.3 ms | 3.4 ms | 3.5 ms | 3.0 ms (reference) |
| ap-gauss-3 | 10 gaussian linear | 26.8 ms | 1.7 ms | 1.8 ms | 6.0 ms (reference) |
| ap-mod-3 | 6 linear, p=2^61-1 | UNSUPPORTED | 414 µs | 409 µs | 195 µs (maxima) |
