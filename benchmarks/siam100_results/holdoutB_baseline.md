# SIAM 100-digit scoring

Cell: correct digits / ms (`~` below need, `✗` wrong, — none).

| experiment | cases | Mathilda on target |
|---|---|---|
| 226-siam100b-c1-oscillatory-tails | 20 | 0 |
| 227-siam100b-c2-billiards | 20 | 0 |
| 228-siam100b-c3-matrix-norms | 20 | 0 |
| 229-siam100b-c4-global-minima | 20 | 5 |
| 230-siam100b-c5-discrete-minimax | 20 | 7 |
| 231-siam100b-c6-lattice-walks | 20 | 2 |
| 232-siam100b-c7-sparse-solves | 20 | 0 |
| 233-siam100b-c8-heat-times | 20 | 0 |
| 234-siam100b-c9-parametric-optima | 20 | 1 |
| 235-siam100b-c10-exit-probabilities | 20 | 0 |
| 236-siam100b-d1-digit-restricted-series | 20 | 16 |
| 237-siam100b-d2-digit-length-series | 20 | 2 |
| 238-siam100b-d3-minimal-energy | 20 | 1 |
| 239-siam100b-d4-series-roots | 20 | 2 |
| 240-siam100b-d5-matrix-polynomial-norms | 20 | 9 |
| 241-siam100b-d6-oscillatory-towers | 20 | 0 |
| 242-siam100b-d7-region-areas | 20 | 0 |
| 243-siam100b-d8-rational-rn-integrals | 20 | 0 |
| 244-siam100b-d9-lattice-walks | 20 | 0 |
| 245-siam100b-d10-random-recurrences | 20 | 0 |
| 246-siam100b-d11-vector-shooting | 20 | 0 |
| 247-siam100b-d12-event-time-root | 20 | 1 |
| 248-siam100b-d13-initial-value-threshold | 20 | 1 |
| 249-siam100b-d14-fem-parameter | 20 | 3 |
| 250-siam100b-d15-blowup-times | 20 | 4 |
| 251-siam100b-d16-heat-functionals | 20 | 8 |
| 252-siam100b-d17-nodal-lines | 20 | 0 |
| 253-siam100b-d18-dirichlet-eigenvalue | 20 | 0 |
| 254-siam100b-d19-design | 20 | 4 |
| 255-siam100b-d20-functional-roots | 20 | 0 |

| system | correct | inexact | wrong | none |
|---|---|---|---|---|
| Mathilda | 0 | 0 | 0 | 600 |
| WL | 275 | 213 | 61 | 51 |
| Py | 318 | 69 | 44 | 169 |
| Maxima | 185 | 26 | 11 | 378 |
| FriCAS | 108 | 16 | 3 | 473 |

### 226-siam100b-c1-oscillatory-tails

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bc1.1 machine | 10 | — | 11 d / 9.5 ms | ~7 d / 5.7 ms | ~7 d / 74 ms | ~5 d / 191 ms | FEWER-DIGITS(0<11) |
| bc1.1 wp50 | 48 | — | ~41 d / 198 ms | ~2 d / 65 ms | — | ~5 d / 42377 ms | FEWER-DIGITS(0<41) |
| bc1.2 machine | 10 | — | 11 d / 11 ms | ✗0 d / 6.4 ms | ✗0 d / 94 ms | ✗0 d / 151 ms | FEWER-DIGITS(0<11) |
| bc1.2 wp50 | 48 | — | ~43 d / 144 ms | ✗0 d / 36 ms | — | ✗0 d / 44153 ms | FEWER-DIGITS(0<43) |
| bc1.3 machine | 10 | — | ~8 d / 11 ms | ~3 d / 8.5 ms | ~3 d / 135 ms | ~4 d / 307 ms | FEWER-DIGITS(0<8) |
| bc1.3 wp50 | 48 | — | ~13 d / 149 ms | ~2 d / 42 ms | — | ~4 d / 90963 ms | FEWER-DIGITS(0<13) |
| bc1.4 machine | 10 | — | ~7 d / 3.7 ms | ✗0 d / 2.4 ms | ✗0 d / 35 ms | — | FEWER-DIGITS(0<7) |
| bc1.4 wp50 | 48 | — | 50 d / 32 ms | 50 d / 1265 ms | — | — | FEWER-DIGITS(0<50) |
| bc1.5 machine | 10 | — | ~7 d / 4.3 ms | — | — | — | FEWER-DIGITS(0<7) |
| bc1.5 wp50 | 48 | — | 50 d / 26 ms | ~2 d / 2602 ms | — | — | FEWER-DIGITS(0<50) |
| bc1.6 machine | 10 | — | ~7 d / 16 ms | ~4 d / 7.2 ms | ~5 d / 121 ms | ~1 d / 453 ms | FEWER-DIGITS(0<7) |
| bc1.6 wp50 | 48 | — | ~7 d / 135 ms | ~1 d / 35 ms | — | ✗0 d / 50087 ms | FEWER-DIGITS(0<7) |
| bc1.7 machine | 10 | — | ~7 d / 4 ms | 10 d / 0.2 ms | ~7 d / 4 ms | — | FEWER-DIGITS(0<10) |
| bc1.7 wp50 | 48 | — | 50 d / 17 ms | 50 d / 548 ms | — | — | FEWER-DIGITS(0<50) |
| bc1.8 machine | 10 | — | 10 d / 291 ms | ✗0 d / 8.3 ms | ✗0 d / 133 ms | — | FEWER-DIGITS(0<10) |
| bc1.8 wp50 | 48 | — | 50 d / 1919 ms | ~8 d / 3816 ms | — | — | FEWER-DIGITS(0<50) |
| bc1.9 machine | 10 | — | ~7 d / 3.5 ms | ✗0 d / 3.8 ms | ✗0 d / 37 ms | — | FEWER-DIGITS(0<7) |
| bc1.9 wp50 | 48 | — | 50 d / 24 ms | 50 d / 506 ms | — | — | FEWER-DIGITS(0<50) |
| bc1.10 machine | 10 | — | ~6 d / 4.3 ms | ~9 d / 0.085 ms | ~6 d / 2.4 ms | — | FEWER-DIGITS(0<9) |
| bc1.10 wp50 | 48 | — | 49 d / 56 ms | 50 d / 3029 ms | — | — | FEWER-DIGITS(0<50) |

### 227-siam100b-c2-billiards

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bc2.1 machine | 10 | — | 12 d / 18 ms | 11 d / 2.3 ms | 11 d / 32 ms | 11 d / 284 ms | FEWER-DIGITS(0<12) |
| bc2.1 wp50 | 48 | — | 50 d / 43 ms | 50 d / 21 ms | 50 d / 90 ms | 50 d / 246 ms | FEWER-DIGITS(0<50) |
| bc2.2 machine | 9 | — | 10 d / 13 ms | 10 d / 1.2 ms | 10 d / 7 ms | 10 d / 14 ms | FEWER-DIGITS(0<10) |
| bc2.2 wp50 | 48 | — | 50 d / 28 ms | 50 d / 13 ms | 50 d / 61 ms | 50 d / 51 ms | FEWER-DIGITS(0<50) |
| bc2.3 machine | 9 | — | 10 d / 22 ms | 10 d / 2 ms | 10 d / 8.8 ms | 10 d / 20 ms | FEWER-DIGITS(0<10) |
| bc2.3 wp50 | 48 | — | 50 d / 43 ms | 50 d / 22 ms | 50 d / 81 ms | 50 d / 56 ms | FEWER-DIGITS(0<50) |
| bc2.4 machine | 10 | — | 12 d / 27 ms | 12 d / 2.7 ms | 12 d / 2.5 ms | 12 d / 22 ms | FEWER-DIGITS(0<12) |
| bc2.4 wp50 | 48 | — | 50 d / 34 ms | 50 d / 6.1 ms | 50 d / 71 ms | 50 d / 54 ms | FEWER-DIGITS(0<50) |
| bc2.5 machine | 10 | — | 11 d / 23 ms | 11 d / 2.4 ms | 11 d / 2.4 ms | 11 d / 20 ms | FEWER-DIGITS(0<11) |
| bc2.5 wp50 | 48 | — | 50 d / 31 ms | 50 d / 6.4 ms | 50 d / 70 ms | 50 d / 51 ms | FEWER-DIGITS(0<50) |
| bc2.6 machine | 10 | — | 10 d / 11 ms | 11 d / 1 ms | 11 d / 2.3 ms | 12 d / 12 ms | FEWER-DIGITS(0<12) |
| bc2.6 wp50 | 48 | — | ~46 d / 17 ms | 50 d / 9.6 ms | 50 d / 47 ms | 50 d / 35 ms | FEWER-DIGITS(0<50) |
| bc2.7 machine | 10 | — | 12 d / 25 ms | 12 d / 2.6 ms | 12 d / 2.3 ms | 12 d / 21 ms | FEWER-DIGITS(0<12) |
| bc2.7 wp50 | 48 | — | 50 d / 32 ms | 50 d / 5.8 ms | 50 d / 74 ms | 50 d / 53 ms | FEWER-DIGITS(0<50) |
| bc2.8 machine | 10 | — | 14 d / 62 ms | 14 d / 6.4 ms | 14 d / 5.6 ms | 14 d / 49 ms | FEWER-DIGITS(0<14) |
| bc2.8 wp50 | 48 | — | ~45 d / 75 ms | 50 d / 9.2 ms | 50 d / 108 ms | 50 d / 103 ms | FEWER-DIGITS(0<50) |
| bc2.9 machine | 10 | — | 11 d / 21 ms | 11 d / 2.1 ms | 11 d / 8.7 ms | 11 d / 22 ms | FEWER-DIGITS(0<11) |
| bc2.9 wp50 | 48 | — | 50 d / 44 ms | 50 d / 23 ms | 50 d / 90 ms | 50 d / 63 ms | FEWER-DIGITS(0<50) |
| bc2.10 machine | 6 | — | 6 d / 49 ms | 6 d / 5.1 ms | 6 d / 4.6 ms | 6 d / 41 ms | FEWER-DIGITS(0<6) |
| bc2.10 wp50 | 47 | — | ~42 d / 63 ms | 50 d / 11 ms | 50 d / 144 ms | 50 d / 103 ms | FEWER-DIGITS(0<50) |

### 228-siam100b-c3-matrix-norms

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bc3.1 machine | 10 | — | 15 d / 236 ms | 15 d / 85 ms | 15 d / 1621 ms | — | FEWER-DIGITS(0<15) |
| bc3.1 wp50 | 48 | — | 50 d / 94 ms | 50 d / 844 ms | — | — | FEWER-DIGITS(0<50) |
| bc3.2 machine | 10 | — | 15 d / 1.9 ms | 15 d / 3.3 ms | 15 d / 11 ms | — | FEWER-DIGITS(0<15) |
| bc3.2 wp50 | 48 | — | 50 d / 32 ms | 50 d / 261 ms | — | — | FEWER-DIGITS(0<50) |
| bc3.3 machine | 10 | — | 15 d / 231 ms | 15 d / 61 ms | 15 d / 1351 ms | — | FEWER-DIGITS(0<15) |
| bc3.3 wp50 | 48 | — | 50 d / 111 ms | 50 d / 991 ms | — | — | FEWER-DIGITS(0<50) |
| bc3.4 machine | 10 | — | 15 d / 211 ms | 15 d / 83 ms | 15 d / 2324 ms | — | FEWER-DIGITS(0<15) |
| bc3.4 wp50 | 48 | — | 50 d / 228 ms | 50 d / 2194 ms | — | — | FEWER-DIGITS(0<50) |
| bc3.5 machine | 10 | — | 15 d / 1286 ms | 15 d / 81 ms | 15 d / 1876 ms | — | FEWER-DIGITS(0<15) |
| bc3.5 wp50 | 48 | — | 50 d / 119 ms | 50 d / 1006 ms | — | — | FEWER-DIGITS(0<50) |
| bc3.6 machine | 10 | — | 15 d / 153 ms | 15 d / 28 ms | 15 d / 850 ms | — | FEWER-DIGITS(0<15) |
| bc3.6 wp50 | 48 | — | 50 d / 53 ms | 50 d / 481 ms | — | — | FEWER-DIGITS(0<50) |
| bc3.7 machine | 10 | — | 15 d / 391 ms | 15 d / 82 ms | — | — | FEWER-DIGITS(0<15) |
| bc3.7 wp50 | 48 | — | 50 d / 219 ms | 50 d / 1118 ms | — | — | FEWER-DIGITS(0<50) |
| bc3.8 machine | 10 | — | 15 d / 106 ms | 15 d / 35 ms | 15 d / 1159 ms | — | FEWER-DIGITS(0<15) |
| bc3.8 wp50 | 48 | — | 50 d / 109 ms | 50 d / 954 ms | — | — | FEWER-DIGITS(0<50) |
| bc3.9 machine | 10 | — | 15 d / 297 ms | 15 d / 8.8 ms | 15 d / 269 ms | — | FEWER-DIGITS(0<15) |
| bc3.9 wp50 | 48 | — | — | 50 d / 12946 ms | — | — | FEWER-DIGITS(0<50) |
| bc3.10 machine | 10 | — | 15 d / 7371 ms | 15 d / 1439 ms | 15 d / 21019 ms | — | FEWER-DIGITS(0<15) |
| bc3.10 wp50 | 48 | — | 50 d / 701 ms | 50 d / 6477 ms | — | — | FEWER-DIGITS(0<50) |

### 229-siam100b-c4-global-minima

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bc4.1 machine | 10 | — | 15 d / 35 ms | 14 d / 41 ms | — | — | FEWER-DIGITS(0<15) |
| bc4.1 wp50 | 48 | — | ✗0 d / 127 ms | 50 d / 46 ms | — | — | FEWER-DIGITS(0<50) |
| bc4.2 machine | 10 | — | ~8 d / 24 ms | ✗0 d / 39 ms | — | — | FEWER-DIGITS(0<8) |
| bc4.2 wp50 | 48 | — | ✗0 d / 43 ms | ✗0 d / 42 ms | — | — | OK(nobody) |
| bc4.3 machine | 10 | — | 15 d / 78 ms | ~1 d / 201 ms | — | — | FEWER-DIGITS(0<15) |
| bc4.3 wp50 | 48 | — | ~1 d / 56 ms | ~1 d / 249 ms | — | — | FEWER-DIGITS(0<1) |
| bc4.4 machine | 10 | — | ~4 d / 901 ms | 14 d / 54 ms | — | — | FEWER-DIGITS(0<14) |
| bc4.4 wp50 | 48 | — | 50 d / 48 ms | 50 d / 61 ms | — | — | FEWER-DIGITS(0<50) |
| bc4.5 machine | 10 | — | ✗0 d / 66 ms | 15 d / 89 ms | — | — | FEWER-DIGITS(0<15) |
| bc4.5 wp50 | 48 | — | ✗0 d / 314 ms | 50 d / 126 ms | — | — | FEWER-DIGITS(0<50) |
| bc4.6 machine | 10 | — | ✗0 d / 28 ms | 13 d / 99 ms | — | — | FEWER-DIGITS(0<13) |
| bc4.6 wp50 | 48 | — | ✗0 d / 139 ms | ~1 d / 104 ms | — | — | FEWER-DIGITS(0<1) |
| bc4.7 machine | 10 | — | ✗0 d / 10 ms | 13 d / 15 ms | — | — | FEWER-DIGITS(0<13) |
| bc4.7 wp50 | 48 | — | 50 d / 162 ms | 50 d / 17 ms | — | — | FEWER-DIGITS(0<50) |
| bc4.8 machine | 10 | — | ✗0 d / 29 ms | ✗0 d / 61 ms | — | — | OK(nobody) |
| bc4.8 wp50 | 48 | — | ✗0 d / 43 ms | ✗0 d / 70 ms | — | — | OK(nobody) |
| bc4.9 machine | 10 | — | ✗0 d / 26 ms | ✗0 d / 65 ms | — | — | OK(nobody) |
| bc4.9 wp50 | 48 | — | ✗0 d / 40 ms | ✗0 d / 71 ms | — | — | OK(nobody) |
| bc4.10 machine | 10 | — | ✗0 d / 22 ms | 15 d / 285 ms | — | — | FEWER-DIGITS(0<15) |
| bc4.10 wp50 | 48 | — | ✗0 d / 54 ms | 50 d / 308 ms | — | — | FEWER-DIGITS(0<50) |

### 230-siam100b-c5-discrete-minimax

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bc5.1 machine | 10 | — | 14 d / 169 ms | 12 d / 321 ms | ✗0 d / 154 ms | — | FEWER-DIGITS(0<14) |
| bc5.1 wp50 | 48 | — | ✗0 d / 1511 ms | — | — | — | OK(nobody) |
| bc5.2 machine | 10 | — | 15 d / 16 ms | 11 d / 1497 ms | ✗0 d / 366 ms | — | FEWER-DIGITS(0<15) |
| bc5.2 wp50 | 48 | — | ✗0 d / 5450 ms | — | — | — | OK(nobody) |
| bc5.3 machine | 10 | — | 15 d / 9.2 ms | 12 d / 1011 ms | ✗0 d / 153 ms | — | FEWER-DIGITS(0<15) |
| bc5.3 wp50 | 48 | — | ~1 d / 2117 ms | — | — | — | FEWER-DIGITS(0<1) |
| bc5.4 machine | 10 | — | 15 d / 26 ms | 11 d / 6020 ms | ✗0 d / 1138 ms | — | FEWER-DIGITS(0<15) |
| bc5.4 wp50 | 48 | — | ✗0 d / 3436 ms | — | — | — | OK(nobody) |
| bc5.5 machine | 10 | — | 14 d / 7.5 ms | 12 d / 1182 ms | ✗0 d / 170 ms | — | FEWER-DIGITS(0<14) |
| bc5.5 wp50 | 48 | — | ✗0 d / 1915 ms | — | — | — | OK(nobody) |
| bc5.6 machine | 10 | — | 10 d / 19 ms | 12 d / 3231 ms | ✗0 d / 658 ms | — | FEWER-DIGITS(0<12) |
| bc5.6 wp50 | 48 | — | ✗0 d / 6618 ms | — | — | — | OK(nobody) |
| bc5.7 machine | 10 | — | ~9 d / 28 ms | 11 d / 116643 ms | ✗0 d / 1953 ms | — | FEWER-DIGITS(0<11) |
| bc5.7 wp50 | 48 | — | ✗0 d / 5930 ms | — | — | — | OK(nobody) |
| bc5.8 machine | 10 | — | ~5 d / 23 ms | 11 d / 4879 ms | — | — | FEWER-DIGITS(0<11) |
| bc5.8 wp50 | 48 | — | ✗0 d / 4650 ms | — | — | — | OK(nobody) |
| bc5.9 machine | 10 | — | ~5 d / 3017 ms | 12 d / 2645 ms | — | — | FEWER-DIGITS(0<12) |
| bc5.9 wp50 | 48 | — | ~3 d / 5380 ms | — | — | — | FEWER-DIGITS(0<3) |
| bc5.10 machine | 10 | — | ~3 d / 4649 ms | 11 d / 35668 ms | — | — | FEWER-DIGITS(0<11) |
| bc5.10 wp50 | 48 | — | ~1 d / 65714 ms | — | — | — | FEWER-DIGITS(0<1) |

### 231-siam100b-c6-lattice-walks

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bc6.1 machine | 10 | — | ~8 d / 363 ms | 12 d / 1108 ms | 12 d / 45095 ms | — | FEWER-DIGITS(0<12) |
| bc6.1 wp50 | 48 | — | ~7 d / 4363 ms | — | — | — | FEWER-DIGITS(0<7) |
| bc6.2 machine | 10 | — | ~6 d / 328 ms | ~5 d / 24680 ms | — | — | FEWER-DIGITS(0<6) |
| bc6.2 wp50 | 48 | — | ~7 d / 7184 ms | — | — | — | FEWER-DIGITS(0<7) |
| bc6.3 machine | 10 | — | ~6 d / 410 ms | 12 d / 631 ms | 11 d / 22953 ms | — | FEWER-DIGITS(0<12) |
| bc6.3 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bc6.4 machine | 10 | — | ~8 d / 281 ms | 12 d / 1513 ms | 11 d / 74255 ms | — | FEWER-DIGITS(0<12) |
| bc6.4 wp50 | 48 | — | ~7 d / 4621 ms | — | — | — | FEWER-DIGITS(0<7) |
| bc6.5 machine | 10 | — | ~7 d / 183 ms | 12 d / 6706 ms | — | — | FEWER-DIGITS(0<12) |
| bc6.5 wp50 | 48 | — | ~6 d / 36880 ms | — | — | — | FEWER-DIGITS(0<6) |
| bc6.6 machine | 10 | — | ~7 d / 226 ms | 11 d / 1453 ms | 11 d / 69125 ms | — | FEWER-DIGITS(0<11) |
| bc6.6 wp50 | 48 | — | ~6 d / 5505 ms | — | — | — | FEWER-DIGITS(0<6) |
| bc6.7 machine | 10 | — | ~5 d / 36 ms | 13 d / 48 ms | 13 d / 1457 ms | — | FEWER-DIGITS(0<13) |
| bc6.7 wp50 | 48 | — | ~8 d / 342 ms | — | — | — | FEWER-DIGITS(0<8) |
| bc6.8 machine | 10 | — | ~6 d / 2171 ms | 14 d / 23176 ms | — | — | FEWER-DIGITS(0<14) |
| bc6.8 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bc6.9 machine | 10 | — | ~7 d / 87 ms | 12 d / 126 ms | 12 d / 3455 ms | — | FEWER-DIGITS(0<12) |
| bc6.9 wp50 | 48 | — | ~7 d / 587 ms | — | — | — | FEWER-DIGITS(0<7) |
| bc6.10 machine | 10 | — | ~7 d / 671 ms | ~7 d / 4136 ms | — | — | FEWER-DIGITS(0<7) |
| bc6.10 wp50 | 48 | — | ~9 d / 8420 ms | — | — | — | FEWER-DIGITS(0<9) |

### 232-siam100b-c7-sparse-solves

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bc7.1 machine | 10 | — | 15 d / 1.7 ms | 15 d / 2.8 ms | 15 d / 14245 ms | 15 d / 1237 ms | FEWER-DIGITS(0<15) |
| bc7.1 wp50 | 48 | — | 50 d / 15 ms | — | — | — | FEWER-DIGITS(0<50) |
| bc7.2 machine | 10 | — | 15 d / 42 ms | 15 d / 21 ms | — | — | FEWER-DIGITS(0<15) |
| bc7.2 wp50 | 48 | — | 50 d / 252 ms | — | — | — | FEWER-DIGITS(0<50) |
| bc7.3 machine | 10 | — | 15 d / 86112 ms | 15 d / 1180 ms | — | — | FEWER-DIGITS(0<15) |
| bc7.3 wp50 | 48 | — | 50 d / 71696 ms | — | — | — | FEWER-DIGITS(0<50) |
| bc7.4 machine | 10 | — | 15 d / 358 ms | 15 d / 743 ms | — | 14 d / 99878 ms | FEWER-DIGITS(0<15) |
| bc7.4 wp50 | 48 | — | 50 d / 488 ms | — | — | — | FEWER-DIGITS(0<50) |
| bc7.5 machine | 10 | — | 15 d / 137 ms | 15 d / 636 ms | — | 14 d / 70525 ms | FEWER-DIGITS(0<15) |
| bc7.5 wp50 | 48 | — | 50 d / 198 ms | — | — | — | FEWER-DIGITS(0<50) |
| bc7.6 machine | 10 | — | 15 d / 1176 ms | 15 d / 13822 ms | — | — | FEWER-DIGITS(0<15) |
| bc7.6 wp50 | 48 | — | 50 d / 1544 ms | — | — | — | FEWER-DIGITS(0<50) |
| bc7.7 machine | 10 | — | 15 d / 21 ms | 15 d / 28 ms | — | 14 d / 45961 ms | FEWER-DIGITS(0<15) |
| bc7.7 wp50 | 48 | — | 50 d / 75 ms | — | — | — | FEWER-DIGITS(0<50) |
| bc7.8 machine | 10 | — | 15 d / 632 ms | 15 d / 5101 ms | — | — | FEWER-DIGITS(0<15) |
| bc7.8 wp50 | 48 | — | 50 d / 1359 ms | — | — | — | FEWER-DIGITS(0<50) |
| bc7.9 machine | 10 | — | 11 d / 106 ms | 12 d / 76 ms | — | — | FEWER-DIGITS(0<12) |
| bc7.9 wp50 | 48 | — | 50 d / 828 ms | — | — | — | FEWER-DIGITS(0<50) |
| bc7.10 machine | 10 | — | 15 d / 3.2 ms | 15 d / 1.5 ms | — | 15 d / 60158 ms | FEWER-DIGITS(0<15) |
| bc7.10 wp50 | 48 | — | 50 d / 44 ms | — | — | — | FEWER-DIGITS(0<50) |

### 233-siam100b-c8-heat-times

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bc8.1 machine | 10 | — | 15 d / 1006 ms | 15 d / 6.9 ms | 15 d / 5.3 ms | — | FEWER-DIGITS(0<15) |
| bc8.1 wp50 | 48 | — | 50 d / 1676 ms | 50 d / 19 ms | 50 d / 168 ms | — | FEWER-DIGITS(0<50) |
| bc8.2 machine | 10 | — | 15 d / 1086 ms | 15 d / 2.7 ms | 15 d / 3.5 ms | — | FEWER-DIGITS(0<15) |
| bc8.2 wp50 | 48 | — | 50 d / 1795 ms | 50 d / 6.2 ms | 50 d / 155 ms | — | FEWER-DIGITS(0<50) |
| bc8.3 machine | 10 | — | 15 d / 976 ms | 15 d / 3.1 ms | 15 d / 4 ms | — | FEWER-DIGITS(0<15) |
| bc8.3 wp50 | 48 | — | 50 d / 1623 ms | 50 d / 5.7 ms | 50 d / 208 ms | — | FEWER-DIGITS(0<50) |
| bc8.4 machine | 10 | — | — | 15 d / 149 ms | — | — | FEWER-DIGITS(0<15) |
| bc8.4 wp50 | 48 | — | — | 50 d / 615 ms | — | — | FEWER-DIGITS(0<50) |
| bc8.5 machine | 10 | — | 15 d / 170 ms | 15 d / 2.2 ms | 15 d / 2.9 ms | — | FEWER-DIGITS(0<15) |
| bc8.5 wp50 | 48 | — | 50 d / 607 ms | 50 d / 14 ms | 50 d / 135 ms | — | FEWER-DIGITS(0<50) |
| bc8.6 machine | 10 | — | 15 d / 83 ms | 15 d / 2.3 ms | 15 d / 3.7 ms | — | FEWER-DIGITS(0<15) |
| bc8.6 wp50 | 48 | — | 50 d / 138 ms | 50 d / 5.1 ms | 50 d / 162 ms | — | FEWER-DIGITS(0<50) |
| bc8.7 machine | 10 | — | 15 d / 1814 ms | 15 d / 4.6 ms | 15 d / 5.3 ms | — | FEWER-DIGITS(0<15) |
| bc8.7 wp50 | 48 | — | 50 d / 2954 ms | 50 d / 22 ms | 50 d / 346 ms | — | FEWER-DIGITS(0<50) |
| bc8.8 machine | 10 | — | 15 d / 3118 ms | 15 d / 9.1 ms | 15 d / 7.1 ms | — | FEWER-DIGITS(0<15) |
| bc8.8 wp50 | 48 | — | 50 d / 5039 ms | 50 d / 143 ms | 50 d / 909 ms | — | FEWER-DIGITS(0<50) |
| bc8.9 machine | 10 | — | 15 d / 19944 ms | 15 d / 4.2 ms | 15 d / 4 ms | — | FEWER-DIGITS(0<15) |
| bc8.9 wp50 | 48 | — | 50 d / 34574 ms | 50 d / 35 ms | 50 d / 222 ms | — | FEWER-DIGITS(0<50) |
| bc8.10 machine | 10 | — | 11 d / 1482 ms | 12 d / 23 ms | 10 d / 6 ms | — | FEWER-DIGITS(0<12) |
| bc8.10 wp50 | 48 | — | 50 d / 2333 ms | 50 d / 23 ms | 50 d / 204 ms | — | FEWER-DIGITS(0<50) |

### 234-siam100b-c9-parametric-optima

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bc9.1 machine | 10 | — | ~7 d / 22567 ms | ~2 d / 5524 ms | — | — | FEWER-DIGITS(0<7) |
| bc9.1 wp50 | 48 | — | ~30 d / 3198 ms | — | — | — | FEWER-DIGITS(0<30) |
| bc9.2 machine | 10 | — | ~6 d / 32259 ms | ~2 d / 11558 ms | — | — | FEWER-DIGITS(0<6) |
| bc9.2 wp50 | 48 | — | ~29 d / 3933 ms | — | — | — | FEWER-DIGITS(0<29) |
| bc9.3 machine | 10 | — | ~7 d / 8734 ms | ✗0 d / 16256 ms | — | — | FEWER-DIGITS(0<7) |
| bc9.3 wp50 | 48 | — | ~30 d / 217 ms | — | — | — | FEWER-DIGITS(0<30) |
| bc9.4 machine | 10 | — | ~9 d / 6003 ms | ~8 d / 30 ms | — | — | FEWER-DIGITS(0<9) |
| bc9.4 wp50 | 48 | — | ~30 d / 333 ms | 50 d / 963 ms | — | — | FEWER-DIGITS(0<50) |
| bc9.5 machine | 10 | — | ~8 d / 10921 ms | ~7 d / 25 ms | — | — | FEWER-DIGITS(0<8) |
| bc9.5 wp50 | 48 | — | ~29 d / 310 ms | 50 d / 3490 ms | — | — | FEWER-DIGITS(0<50) |
| bc9.6 machine | 10 | — | 10 d / 12874 ms | ~8 d / 35 ms | — | — | FEWER-DIGITS(0<10) |
| bc9.6 wp50 | 48 | — | ~30 d / 651 ms | — | — | — | FEWER-DIGITS(0<30) |
| bc9.7 machine | 10 | — | ~9 d / 8246 ms | ~7 d / 44 ms | — | — | FEWER-DIGITS(0<9) |
| bc9.7 wp50 | 48 | — | ~30 d / 643 ms | 50 d / 5087 ms | — | — | FEWER-DIGITS(0<50) |
| bc9.8 machine | 10 | — | 10 d / 8905 ms | ~7 d / 18 ms | — | — | FEWER-DIGITS(0<10) |
| bc9.8 wp50 | 48 | — | ~30 d / 629 ms | 50 d / 220 ms | — | — | FEWER-DIGITS(0<50) |
| bc9.9 machine | 10 | — | ~4 d / 42395 ms | — | — | — | FEWER-DIGITS(0<4) |
| bc9.9 wp50 | 48 | — | ~30 d / 12593 ms | — | — | — | FEWER-DIGITS(0<30) |
| bc9.10 machine | 10 | — | — | ✗0 d / 879 ms | — | — | OK(nobody) |
| bc9.10 wp50 | 48 | — | ~4 d / 3928 ms | — | — | — | FEWER-DIGITS(0<4) |

### 235-siam100b-c10-exit-probabilities

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bc10.1 machine | 10 | — | 15 d / 440 ms | 15 d / 0.65 ms | 15 d / 0.35 ms | 15 d / 22 ms | FEWER-DIGITS(0<15) |
| bc10.1 wp50 | 48 | — | 50 d / 53 ms | 50 d / 1.8 ms | 50 d / 1.6 ms | 50 d / 14 ms | FEWER-DIGITS(0<50) |
| bc10.2 machine | 10 | — | 15 d / 555 ms | 15 d / 0.42 ms | 15 d / 0.49 ms | 15 d / 9.5 ms | FEWER-DIGITS(0<15) |
| bc10.2 wp50 | 48 | — | 50 d / 86 ms | 50 d / 1.9 ms | 50 d / 3.1 ms | 50 d / 11 ms | FEWER-DIGITS(0<50) |
| bc10.3 machine | 10 | — | 15 d / 342 ms | 15 d / 0.43 ms | 15 d / 0.67 ms | 15 d / 8.5 ms | FEWER-DIGITS(0<15) |
| bc10.3 wp50 | 48 | — | 50 d / 97 ms | 50 d / 1.9 ms | 50 d / 4.3 ms | 50 d / 11 ms | FEWER-DIGITS(0<50) |
| bc10.4 machine | 10 | — | — | 15 d / 4.5 ms | 15 d / 140 ms | 15 d / 15 ms | FEWER-DIGITS(0<15) |
| bc10.4 wp50 | 48 | — | — | 50 d / 50 ms | 50 d / 280 ms | 50 d / 60 ms | FEWER-DIGITS(0<50) |
| bc10.5 machine | 10 | — | 15 d / 8307 ms | 15 d / 0.48 ms | 15 d / 0.69 ms | 15 d / 8.1 ms | FEWER-DIGITS(0<15) |
| bc10.5 wp50 | 48 | — | 50 d / 108 ms | 50 d / 2.4 ms | 50 d / 1.8 ms | 50 d / 11 ms | FEWER-DIGITS(0<50) |
| bc10.6 machine | 10 | — | 15 d / 2195 ms | 15 d / 0.57 ms | 15 d / 0.29 ms | 15 d / 7 ms | FEWER-DIGITS(0<15) |
| bc10.6 wp50 | 48 | — | 50 d / 1483 ms | 50 d / 0.62 ms | 50 d / 1.4 ms | 50 d / 10 ms | FEWER-DIGITS(0<50) |
| bc10.7 machine | 10 | — | 15 d / 179 ms | 15 d / 0.36 ms | 15 d / 0.26 ms | 15 d / 4.9 ms | FEWER-DIGITS(0<15) |
| bc10.7 wp50 | 48 | — | 50 d / 51 ms | 50 d / 2.1 ms | 50 d / 0.25 ms | 50 d / 6.3 ms | FEWER-DIGITS(0<50) |
| bc10.8 machine | 10 | — | — | 15 d / 19 ms | — | — | FEWER-DIGITS(0<15) |
| bc10.8 wp50 | 48 | — | — | 50 d / 41 ms | — | — | FEWER-DIGITS(0<50) |
| bc10.9 machine | 10 | — | 15 d / 3548 ms | 15 d / 2.7 ms | ~5 d / 5.4 ms | ~5 d / 7.9 ms | FEWER-DIGITS(0<15) |
| bc10.9 wp50 | 48 | — | ~15 d / 109 ms | 50 d / 25 ms | 50 d / 1403 ms | 50 d / 126 ms | FEWER-DIGITS(0<50) |
| bc10.10 machine | 10 | — | — | 15 d / 4.1 ms | 15 d / 4.2 ms | 15 d / 12 ms | FEWER-DIGITS(0<15) |
| bc10.10 wp50 | 48 | — | 50 d / 109487 ms | 50 d / 5.1 ms | 50 d / 8.8 ms | 50 d / 18 ms | FEWER-DIGITS(0<50) |

### 236-siam100b-d1-digit-restricted-series

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd1.1 machine | 10 | — | ✗0 d / 233 ms | ✗0 d / 40 ms | — | — | OK(nobody) |
| bd1.1 wp50 | 48 | — | ✗0 d / 42 ms | ✗0 d / 899 ms | — | — | OK(nobody) |
| bd1.2 machine | 10 | — | ✗0 d / 191 ms | ✗0 d / 0.14 ms | — | — | OK(nobody) |
| bd1.2 wp50 | 48 | — | — | ✗0 d / 0.11 ms | — | — | OK(nobody) |
| bd1.3 machine | 10 | — | ✗0 d / 150 ms | ✗0 d / 0.12 ms | — | — | OK(nobody) |
| bd1.3 wp50 | 48 | — | ✗0 d / 42 ms | ✗0 d / 0.11 ms | — | — | OK(nobody) |
| bd1.4 machine | 10 | — | ✗0 d / 154 ms | ✗0 d / 40 ms | — | — | OK(nobody) |
| bd1.4 wp50 | 48 | — | — | ✗0 d / 852 ms | — | — | OK(nobody) |
| bd1.5 machine | 10 | — | ✗0 d / 151 ms | ✗0 d / 0.11 ms | — | — | OK(nobody) |
| bd1.5 wp50 | 48 | — | ✗0 d / 43 ms | ✗0 d / 0.081 ms | — | — | OK(nobody) |
| bd1.6 machine | 10 | — | ~3 d / 202 ms | ~2 d / 46 ms | — | — | FEWER-DIGITS(0<3) |
| bd1.6 wp50 | 48 | — | — | ~2 d / 42 ms | — | — | FEWER-DIGITS(0<2) |
| bd1.7 machine | 10 | — | ✗0 d / 161 ms | ✗0 d / 0.35 ms | — | — | OK(nobody) |
| bd1.7 wp50 | 48 | — | — | ✗0 d / 0.34 ms | — | — | OK(nobody) |
| bd1.8 machine | 10 | — | ✗0 d / 182 ms | ✗0 d / 0.33 ms | — | — | OK(nobody) |
| bd1.8 wp50 | 48 | — | ✗0 d / 43 ms | ✗0 d / 0.34 ms | — | — | OK(nobody) |
| bd1.9 machine | 10 | — | ~9 d / 148 ms | ~2 d / 0.35 ms | — | — | FEWER-DIGITS(0<9) |
| bd1.9 wp50 | 48 | — | ~9 d / 42 ms | ~2 d / 0.35 ms | — | — | FEWER-DIGITS(0<9) |
| bd1.10 machine | 10 | — | ✗0 d / 157 ms | ✗0 d / 6.7 ms | — | — | OK(nobody) |
| bd1.10 wp50 | 48 | — | — | ✗0 d / 7.2 ms | — | — | OK(nobody) |

### 237-siam100b-d2-digit-length-series

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd2.1 machine | 10 | — | ~3 d / 310 ms | ✗0 d / 39 ms | — | — | FEWER-DIGITS(0<3) |
| bd2.1 wp50 | 48 | — | ~4 d / 151 ms | ~1 d / 876 ms | — | — | FEWER-DIGITS(0<4) |
| bd2.2 machine | 10 | — | ~5 d / 165 ms | ~1 d / 40 ms | — | — | FEWER-DIGITS(0<5) |
| bd2.2 wp50 | 48 | — | ~4 d / 133 ms | ~1 d / 839 ms | — | — | FEWER-DIGITS(0<4) |
| bd2.3 machine | 10 | — | 10 d / 16 ms | ~1 d / 41 ms | — | — | FEWER-DIGITS(0<10) |
| bd2.3 wp50 | 48 | — | ~26 d / 60 ms | ~1 d / 860 ms | — | — | FEWER-DIGITS(0<26) |
| bd2.4 machine | 10 | — | ~2 d / 25 ms | ~1 d / 42 ms | — | — | FEWER-DIGITS(0<2) |
| bd2.4 wp50 | 48 | — | ~2 d / 231 ms | ~1 d / 877 ms | — | — | FEWER-DIGITS(0<2) |
| bd2.5 machine | 10 | — | ✗0 d / 158 ms | ✗0 d / 40 ms | — | — | OK(nobody) |
| bd2.5 wp50 | 48 | — | ✗0 d / 81 ms | ✗0 d / 825 ms | — | — | OK(nobody) |
| bd2.6 machine | 10 | — | ~4 d / 162 ms | ✗0 d / 39 ms | — | — | FEWER-DIGITS(0<4) |
| bd2.6 wp50 | 48 | — | ~4 d / 106 ms | ~1 d / 845 ms | — | — | FEWER-DIGITS(0<4) |
| bd2.7 machine | 10 | — | ~1 d / 158 ms | ✗0 d / 40 ms | — | — | FEWER-DIGITS(0<1) |
| bd2.7 wp50 | 48 | — | ~1 d / 71 ms | ✗0 d / 848 ms | — | — | FEWER-DIGITS(0<1) |
| bd2.8 machine | 10 | — | ~2 d / 154 ms | ~1 d / 41 ms | — | — | FEWER-DIGITS(0<2) |
| bd2.8 wp50 | 48 | — | ~3 d / 70 ms | ~1 d / 837 ms | — | — | FEWER-DIGITS(0<3) |
| bd2.9 machine | 10 | — | ✗0 d / 154 ms | ~2 d / 39 ms | — | — | FEWER-DIGITS(0<2) |
| bd2.9 wp50 | 48 | — | ✗0 d / 73 ms | ~1 d / 829 ms | — | — | FEWER-DIGITS(0<1) |
| bd2.10 machine | 10 | — | ~1 d / 158 ms | ✗0 d / 40 ms | — | — | FEWER-DIGITS(0<1) |
| bd2.10 wp50 | 48 | — | ~1 d / 77 ms | ✗0 d / 860 ms | — | — | FEWER-DIGITS(0<1) |

### 238-siam100b-d3-minimal-energy

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd3.1 machine | 10 | — | 15 d / 6 ms | 11 d / 17 ms | 15 d / 124 ms | — | FEWER-DIGITS(0<15) |
| bd3.1 wp50 | 48 | — | 50 d / 44 ms | — | — | — | FEWER-DIGITS(0<50) |
| bd3.2 machine | 10 | — | 15 d / 62 ms | 12 d / 48 ms | 15 d / 1002 ms | — | FEWER-DIGITS(0<15) |
| bd3.2 wp50 | 48 | — | 50 d / 209 ms | — | — | — | FEWER-DIGITS(0<50) |
| bd3.3 machine | 10 | — | ✗0 d / 129 ms | 12 d / 452 ms | 15 d / 3621 ms | — | FEWER-DIGITS(0<15) |
| bd3.3 wp50 | 48 | — | ✗0 d / 1011 ms | — | — | — | OK(nobody) |
| bd3.4 machine | 10 | — | 15 d / 103 ms | 12 d / 192 ms | 12 d / 5739 ms | — | FEWER-DIGITS(0<15) |
| bd3.4 wp50 | 48 | — | ~18 d / 1114 ms | — | — | — | FEWER-DIGITS(0<18) |
| bd3.5 machine | 10 | — | 10 d / 189 ms | 11 d / 413 ms | 13 d / 13795 ms | — | FEWER-DIGITS(0<13) |
| bd3.5 wp50 | 48 | — | ~10 d / 2022 ms | — | — | — | FEWER-DIGITS(0<10) |
| bd3.6 machine | 10 | — | 13 d / 78 ms | 11 d / 160 ms | 12 d / 3648 ms | — | FEWER-DIGITS(0<13) |
| bd3.6 wp50 | 48 | — | ~13 d / 690 ms | — | — | — | FEWER-DIGITS(0<13) |
| bd3.7 machine | 10 | — | 15 d / 70 ms | 12 d / 180 ms | 15 d / 3421 ms | — | FEWER-DIGITS(0<15) |
| bd3.7 wp50 | 48 | — | ~16 d / 648 ms | — | — | — | FEWER-DIGITS(0<16) |
| bd3.8 machine | 10 | — | 15 d / 41 ms | 12 d / 244 ms | 15 d / 716 ms | — | FEWER-DIGITS(0<15) |
| bd3.8 wp50 | 48 | — | 50 d / 452 ms | — | — | — | FEWER-DIGITS(0<50) |
| bd3.9 machine | 10 | — | ~6 d / 706 ms | 12 d / 2968 ms | — | — | FEWER-DIGITS(0<12) |
| bd3.9 wp50 | 48 | — | ~6 d / 6969 ms | — | — | — | FEWER-DIGITS(0<6) |
| bd3.10 machine | 10 | — | ~5 d / 200 ms | ~7 d / 689 ms | ~6 d / 12054 ms | — | FEWER-DIGITS(0<7) |
| bd3.10 wp50 | 48 | — | ~5 d / 1526 ms | — | — | — | FEWER-DIGITS(0<5) |

### 239-siam100b-d4-series-roots

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd4.1 machine | 10 | — | 15 d / 2.4 ms | 15 d / 0.41 ms | 15 d / 2.5 ms | — | FEWER-DIGITS(0<15) |
| bd4.1 wp50 | 48 | — | 50 d / 1.4 ms | 50 d / 1.4 ms | 50 d / 77 ms | — | FEWER-DIGITS(0<50) |
| bd4.2 machine | 10 | — | 10 d / 113 ms | — | — | — | FEWER-DIGITS(0<10) |
| bd4.2 wp50 | 48 | — | ~26 d / 412 ms | — | — | — | FEWER-DIGITS(0<26) |
| bd4.3 machine | 10 | — | 10 d / 146 ms | — | — | — | FEWER-DIGITS(0<10) |
| bd4.3 wp50 | 48 | — | ~27 d / 889 ms | — | — | — | FEWER-DIGITS(0<27) |
| bd4.4 machine | 10 | — | 15 d / 61 ms | 15 d / 17 ms | — | — | FEWER-DIGITS(0<15) |
| bd4.4 wp50 | 48 | — | 50 d / 89 ms | 50 d / 31 ms | — | — | FEWER-DIGITS(0<50) |
| bd4.5 machine | 10 | — | — | ✗0 d / 72 ms | — | — | OK(nobody) |
| bd4.5 wp50 | 48 | — | — | ✗0 d / 105 ms | — | — | OK(nobody) |
| bd4.6 machine | 10 | — | ~6 d / 3447 ms | — | — | — | FEWER-DIGITS(0<6) |
| bd4.6 wp50 | 48 | — | ~6 d / 4410 ms | — | — | — | FEWER-DIGITS(0<6) |
| bd4.7 machine | 10 | — | 15 d / 31 ms | 15 d / 16 ms | — | — | FEWER-DIGITS(0<15) |
| bd4.7 wp50 | 48 | — | 50 d / 61 ms | 50 d / 17 ms | — | — | FEWER-DIGITS(0<50) |
| bd4.8 machine | 10 | — | 15 d / 3.6 ms | 15 d / 71 ms | — | — | FEWER-DIGITS(0<15) |
| bd4.8 wp50 | 48 | — | 50 d / 61 ms | 50 d / 791 ms | — | — | FEWER-DIGITS(0<50) |
| bd4.9 machine | 10 | — | — | 15 d / 6220 ms | — | — | FEWER-DIGITS(0<15) |
| bd4.9 wp50 | 48 | — | 50 d / 4044 ms | 50 d / 7278 ms | — | — | FEWER-DIGITS(0<50) |
| bd4.10 machine | 10 | — | 10 d / 231 ms | — | — | — | FEWER-DIGITS(0<10) |
| bd4.10 wp50 | 48 | — | ~28 d / 748 ms | — | — | — | FEWER-DIGITS(0<28) |

### 240-siam100b-d5-matrix-polynomial-norms

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd5.1 machine | 10 | — | ✗0 d / 64 ms | ✗0 d / 28 ms | — | — | OK(nobody) |
| bd5.1 wp50 | 48 | — | ✗0 d / 243 ms | — | — | — | OK(nobody) |
| bd5.2 machine | 10 | — | 13 d / 1819 ms | 10 d / 36 ms | 15 d / 398 ms | — | FEWER-DIGITS(0<15) |
| bd5.2 wp50 | 48 | — | 50 d / 5736 ms | — | — | — | FEWER-DIGITS(0<50) |
| bd5.3 machine | 10 | — | 10 d / 651 ms | 12 d / 69 ms | — | — | FEWER-DIGITS(0<12) |
| bd5.3 wp50 | 48 | — | 50 d / 50190 ms | — | — | — | FEWER-DIGITS(0<50) |
| bd5.4 machine | 10 | — | ✗0 d / 27 ms | ✗0 d / 33 ms | — | — | OK(nobody) |
| bd5.4 wp50 | 48 | — | ✗0 d / 1325 ms | — | — | — | OK(nobody) |
| bd5.5 machine | 10 | — | ✗0 d / 127 ms | ✗0 d / 64 ms | — | — | OK(nobody) |
| bd5.5 wp50 | 48 | — | ✗0 d / 15669 ms | — | — | — | OK(nobody) |
| bd5.6 machine | 10 | — | 14 d / 104 ms | 14 d / 14 ms | — | — | FEWER-DIGITS(0<14) |
| bd5.6 wp50 | 48 | — | 50 d / 5786 ms | — | — | — | FEWER-DIGITS(0<50) |
| bd5.7 machine | 10 | — | 15 d / 3.2 ms | 12 d / 2.3 ms | — | — | FEWER-DIGITS(0<15) |
| bd5.7 wp50 | 48 | — | 50 d / 1344 ms | — | — | — | FEWER-DIGITS(0<50) |
| bd5.8 machine | 10 | — | ✗0 d / 290 ms | ~3 d / 114 ms | — | — | FEWER-DIGITS(0<3) |
| bd5.8 wp50 | 48 | — | ✗0 d / 24765 ms | — | — | — | OK(nobody) |
| bd5.9 machine | 10 | — | ✗0 d / 34 ms | ✗0 d / 27 ms | — | — | OK(nobody) |
| bd5.9 wp50 | 48 | — | ✗0 d / 2158 ms | — | — | — | OK(nobody) |
| bd5.10 machine | 10 | — | ✗0 d / 107 ms | 13 d / 17 ms | — | — | FEWER-DIGITS(0<13) |
| bd5.10 wp50 | 48 | — | 50 d / 26821 ms | — | — | — | FEWER-DIGITS(0<50) |

### 241-siam100b-d6-oscillatory-towers

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd6.1 machine | 10 | — | 14 d / 2.1 ms | 15 d / 0.009 ms | 15 d / 3.3 ms | 15 d / 7.7 ms | FEWER-DIGITS(0<15) |
| bd6.1 wp50 | 48 | — | 50 d / 21 ms | 50 d / 11 ms | ~7 d / 6.6 ms | 50 d / 1108 ms | FEWER-DIGITS(0<50) |
| bd6.2 machine | 10 | — | 11 d / 5.5 ms | 15 d / 0.17 ms | 15 d / 5.5 ms | 14 d / 28 ms | FEWER-DIGITS(0<15) |
| bd6.2 wp50 | 48 | — | 50 d / 108 ms | 50 d / 42 ms | ~7 d / 284 ms | ~41 d / 8524 ms | FEWER-DIGITS(0<50) |
| bd6.3 machine | 10 | — | ~9 d / 2.3 ms | 15 d / 0.066 ms | 15 d / 2.7 ms | 15 d / 3.2 ms | FEWER-DIGITS(0<15) |
| bd6.3 wp50 | 48 | — | 50 d / 43 ms | 50 d / 32 ms | ~8 d / 32 ms | 50 d / 8949 ms | FEWER-DIGITS(0<50) |
| bd6.4 machine | 10 | — | 10 d / 2.7 ms | 15 d / 0.09 ms | 15 d / 6 ms | 15 d / 4.4 ms | FEWER-DIGITS(0<15) |
| bd6.4 wp50 | 48 | — | 50 d / 71 ms | 50 d / 86 ms | ~6 d / 49 ms | 50 d / 12106 ms | FEWER-DIGITS(0<50) |
| bd6.5 machine | 10 | — | 11 d / 9.2 ms | 15 d / 0.36 ms | 15 d / 7.5 ms | 14 d / 76 ms | FEWER-DIGITS(0<15) |
| bd6.5 wp50 | 48 | — | 50 d / 153 ms | ~10 d / 61 ms | ~6 d / 168 ms | ~31 d / 5919 ms | FEWER-DIGITS(0<50) |
| bd6.6 machine | 10 | — | 13 d / 5.2 ms | 14 d / 1.6 ms | 13 d / 30 ms | 14 d / 104 ms | FEWER-DIGITS(0<14) |
| bd6.6 wp50 | 48 | — | 50 d / 546 ms | ✗0 d / 56 ms | — | ~14 d / 8248 ms | FEWER-DIGITS(0<50) |
| bd6.7 machine | 10 | — | 15 d / 5.1 ms | 14 d / 0.15 ms | 14 d / 4.3 ms | 12 d / 4.9 ms | FEWER-DIGITS(0<15) |
| bd6.7 wp50 | 48 | — | 50 d / 66 ms | 50 d / 28 ms | ~5 d / 18 ms | 50 d / 5240 ms | FEWER-DIGITS(0<50) |
| bd6.8 machine | 10 | — | 13 d / 5 ms | 14 d / 0.19 ms | 14 d / 5.6 ms | 14 d / 16 ms | FEWER-DIGITS(0<14) |
| bd6.8 wp50 | 48 | — | 50 d / 47 ms | 50 d / 40 ms | ~9 d / 360 ms | ~46 d / 11402 ms | FEWER-DIGITS(0<50) |
| bd6.9 machine | 10 | — | 12 d / 5 ms | 15 d / 0.23 ms | 15 d / 4.7 ms | 13 d / 7.6 ms | FEWER-DIGITS(0<15) |
| bd6.9 wp50 | 48 | — | 50 d / 101 ms | 50 d / 68 ms | — | ~31 d / 11118 ms | FEWER-DIGITS(0<50) |
| bd6.10 machine | 10 | — | 11 d / 15 ms | 15 d / 0.37 ms | 15 d / 7.6 ms | 15 d / 52 ms | FEWER-DIGITS(0<15) |
| bd6.10 wp50 | 48 | — | 50 d / 189 ms | ~18 d / 76 ms | ~7 d / 344 ms | ~30 d / 11041 ms | FEWER-DIGITS(0<50) |

### 242-siam100b-d7-region-areas

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd7.1 machine | 10 | — | 15 d / 3.4 ms | 15 d / 0.16 ms | 15 d / 4.6 ms | — | FEWER-DIGITS(0<15) |
| bd7.1 wp50 | 48 | — | 50 d / 7.5 ms | 50 d / 9.4 ms | ~6 d / 6.4 ms | — | FEWER-DIGITS(0<50) |
| bd7.2 machine | 10 | — | 15 d / 3 ms | 15 d / 0.27 ms | 15 d / 5.1 ms | — | FEWER-DIGITS(0<15) |
| bd7.2 wp50 | 48 | — | 50 d / 7.4 ms | 50 d / 4.6 ms | ~9 d / 9.9 ms | — | FEWER-DIGITS(0<50) |
| bd7.3 machine | 10 | — | 15 d / 3 ms | 15 d / 0.45 ms | 15 d / 15 ms | ~9 d / 270 ms | FEWER-DIGITS(0<15) |
| bd7.3 wp50 | 48 | — | 50 d / 82 ms | 50 d / 12 ms | ~4 d / 267 ms | ~7 d / 1973 ms | FEWER-DIGITS(0<50) |
| bd7.4 machine | 10 | — | 14 d / 7.6 ms | 15 d / 1.6 ms | 15 d / 17 ms | — | FEWER-DIGITS(0<15) |
| bd7.4 wp50 | 48 | — | 50 d / 84 ms | 50 d / 28 ms | ~5 d / 91 ms | — | FEWER-DIGITS(0<50) |
| bd7.5 machine | 10 | — | 13 d / 3.3 ms | 15 d / 0.5 ms | 15 d / 6.3 ms | ~9 d / 351 ms | FEWER-DIGITS(0<15) |
| bd7.5 wp50 | 48 | — | 50 d / 24 ms | 50 d / 3.6 ms | ~4 d / 23 ms | ~8 d / 223 ms | FEWER-DIGITS(0<50) |
| bd7.6 machine | 10 | — | 14 d / 2.6 ms | 15 d / 0.25 ms | 15 d / 4.5 ms | — | FEWER-DIGITS(0<15) |
| bd7.6 wp50 | 48 | — | 50 d / 31 ms | 50 d / 11 ms | ~4 d / 61 ms | — | FEWER-DIGITS(0<50) |
| bd7.7 machine | 10 | — | 13 d / 3.7 ms | 15 d / 1.2 ms | 15 d / 7.9 ms | — | FEWER-DIGITS(0<15) |
| bd7.7 wp50 | 48 | — | 50 d / 36 ms | 50 d / 7.9 ms | ~4 d / 134 ms | — | FEWER-DIGITS(0<50) |
| bd7.8 machine | 10 | — | 11 d / 9.5 ms | 12 d / 10 ms | 12 d / 149 ms | — | FEWER-DIGITS(0<12) |
| bd7.8 wp50 | 48 | — | 50 d / 207 ms | ~45 d / 139 ms | ~4 d / 1401 ms | — | FEWER-DIGITS(0<50) |
| bd7.9 machine | 10 | — | 14 d / 6 ms | 15 d / 0.33 ms | 15 d / 9.6 ms | — | FEWER-DIGITS(0<15) |
| bd7.9 wp50 | 48 | — | 50 d / 30 ms | 50 d / 36 ms | ~5 d / 36 ms | — | FEWER-DIGITS(0<50) |
| bd7.10 machine | 10 | — | 13 d / 6846 ms | 14 d / 7 ms | 15 d / 121 ms | — | FEWER-DIGITS(0<15) |
| bd7.10 wp50 | 48 | — | 50 d / 7571 ms | 50 d / 173 ms | ~5 d / 3315 ms | — | FEWER-DIGITS(0<50) |

### 243-siam100b-d8-rational-rn-integrals

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd8.1 machine | 10 | — | 10 d / 2.3 ms | 15 d / 0.52 ms | 15 d / 2.4 ms | — | FEWER-DIGITS(0<15) |
| bd8.1 wp50 | 48 | — | 50 d / 12 ms | 50 d / 59 ms | — | — | FEWER-DIGITS(0<50) |
| bd8.2 machine | 10 | — | 13 d / 2.5 ms | 15 d / 2 ms | 15 d / 3.1 ms | — | FEWER-DIGITS(0<15) |
| bd8.2 wp50 | 48 | — | 50 d / 28 ms | 50 d / 204 ms | — | — | FEWER-DIGITS(0<50) |
| bd8.3 machine | 10 | — | ~5 d / 15 ms | 15 d / 858 ms | 15 d / 947 ms | — | FEWER-DIGITS(0<15) |
| bd8.3 wp50 | 48 | — | ~17 d / 53242 ms | — | — | — | FEWER-DIGITS(0<17) |
| bd8.4 machine | 10 | — | ~8 d / 13 ms | 15 d / 2400 ms | 15 d / 1657 ms | — | FEWER-DIGITS(0<15) |
| bd8.4 wp50 | 48 | — | ~17 d / 60535 ms | — | — | — | FEWER-DIGITS(0<17) |
| bd8.5 machine | 10 | — | ~9 d / 662 ms | — | — | — | FEWER-DIGITS(0<9) |
| bd8.5 wp50 | 48 | — | ~11 d / 42931 ms | — | — | — | FEWER-DIGITS(0<11) |
| bd8.6 machine | 10 | — | ~8 d / 2304 ms | — | — | — | FEWER-DIGITS(0<8) |
| bd8.6 wp50 | 48 | — | ~9 d / 32901 ms | — | — | — | FEWER-DIGITS(0<9) |
| bd8.7 machine | 10 | — | ~7 d / 27 ms | 15 d / 33 ms | 15 d / 718 ms | — | FEWER-DIGITS(0<15) |
| bd8.7 wp50 | 48 | — | ~16 d / 27822 ms | ~40 d / 39380 ms | — | — | FEWER-DIGITS(0<40) |
| bd8.8 machine | 10 | — | ~7 d / 449 ms | 15 d / 10553 ms | — | — | FEWER-DIGITS(0<15) |
| bd8.8 wp50 | 48 | — | ~10 d / 21209 ms | — | — | — | FEWER-DIGITS(0<10) |
| bd8.9 machine | 10 | — | ~5 d / 4857 ms | — | — | — | FEWER-DIGITS(0<5) |
| bd8.9 wp50 | 48 | — | ~5 d / 24358 ms | — | — | — | FEWER-DIGITS(0<5) |
| bd8.10 machine | 10 | — | ~8 d / 24 ms | 15 d / 35 ms | 15 d / 827 ms | — | FEWER-DIGITS(0<15) |
| bd8.10 wp50 | 48 | — | ~16 d / 34170 ms | ~40 d / 60611 ms | — | — | FEWER-DIGITS(0<40) |

### 244-siam100b-d9-lattice-walks

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd9.1 machine | 10 | — | 15 d / 0.073 ms | 15 d / 0.035 ms | 15 d / 15 ms | 15 d / 3.8 ms | FEWER-DIGITS(0<15) |
| bd9.1 wp50 | 48 | — | 50 d / 1.8 ms | 50 d / 13 ms | 50 d / 1.8 ms | 50 d / 0.95 ms | FEWER-DIGITS(0<50) |
| bd9.2 machine | 10 | — | 15 d / 0.5 ms | 15 d / 0.33 ms | 15 d / 1231 ms | 15 d / 207 ms | FEWER-DIGITS(0<15) |
| bd9.2 wp50 | 48 | — | 50 d / 19 ms | 50 d / 184 ms | 50 d / 30 ms | 50 d / 17 ms | FEWER-DIGITS(0<50) |
| bd9.3 machine | 10 | — | 14 d / 1.6 ms | 14 d / 1.3 ms | 14 d / 14720 ms | 15 d / 1677 ms | FEWER-DIGITS(0<15) |
| bd9.3 wp50 | 48 | — | 50 d / 160 ms | 50 d / 1049 ms | 50 d / 264 ms | 50 d / 149 ms | FEWER-DIGITS(0<50) |
| bd9.4 machine | 10 | — | 15 d / 1.6 ms | 15 d / 1 ms | 15 d / 9476 ms | 15 d / 1357 ms | FEWER-DIGITS(0<15) |
| bd9.4 wp50 | 48 | — | 50 d / 67 ms | 50 d / 667 ms | 50 d / 104 ms | 50 d / 59 ms | FEWER-DIGITS(0<50) |
| bd9.5 machine | 10 | — | 15 d / 0.22 ms | 15 d / 0.084 ms | 15 d / 202 ms | 15 d / 37 ms | FEWER-DIGITS(0<15) |
| bd9.5 wp50 | 48 | — | 50 d / 6.4 ms | 50 d / 53 ms | 50 d / 8.7 ms | 50 d / 5.2 ms | FEWER-DIGITS(0<50) |
| bd9.6 machine | 10 | — | 15 d / 0.24 ms | 15 d / 0.12 ms | 15 d / 306 ms | 15 d / 52 ms | FEWER-DIGITS(0<15) |
| bd9.6 wp50 | 48 | — | 50 d / 31 ms | 50 d / 77 ms | 50 d / 44 ms | 50 d / 33 ms | FEWER-DIGITS(0<50) |
| bd9.7 machine | 10 | — | 14 d / 72 ms | 14 d / 48 ms | — | — | FEWER-DIGITS(0<14) |
| bd9.7 wp50 | 48 | — | 50 d / 2690 ms | 50 d / 25109 ms | 50 d / 5271 ms | 50 d / 2184 ms | FEWER-DIGITS(0<50) |
| bd9.8 machine | 10 | — | 15 d / 0.026 ms | 15 d / 0.036 ms | 15 d / 12 ms | 15 d / 1 ms | FEWER-DIGITS(0<15) |
| bd9.8 wp50 | 48 | — | 50 d / 1.4 ms | 50 d / 13 ms | 50 d / 3 ms | 50 d / 1.7 ms | FEWER-DIGITS(0<50) |
| bd9.9 machine | 10 | — | 14 d / 0.063 ms | 15 d / 0.056 ms | 14 d / 266 ms | 14 d / 6.6 ms | FEWER-DIGITS(0<15) |
| bd9.9 wp50 | 48 | — | 50 d / 17 ms | 50 d / 124 ms | 50 d / 42 ms | 50 d / 32 ms | FEWER-DIGITS(0<50) |
| bd9.10 machine | 10 | — | 14 d / 47 ms | 14 d / 36 ms | — | — | FEWER-DIGITS(0<14) |
| bd9.10 wp50 | 48 | — | 50 d / 4442 ms | 50 d / 39746 ms | 50 d / 10700 ms | 50 d / 4751 ms | FEWER-DIGITS(0<50) |

### 245-siam100b-d10-random-recurrences

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd10.1 machine | 2 | — | 3 d / 411 ms | 3 d / 178 ms | 3 d / 6349 ms | 3 d / 146 ms | FEWER-DIGITS(0<3) |
| bd10.1 wp50 | 1 | — | 2 d / 28 ms | 2 d / 5.5 ms | 2 d / 74 ms | 2 d / 18 ms | FEWER-DIGITS(0<2) |
| bd10.2 machine | 3 | — | 5 d / 406 ms | 4 d / 179 ms | 5 d / 6372 ms | 3 d / 139 ms | FEWER-DIGITS(0<5) |
| bd10.2 wp50 | 2 | — | 3 d / 38 ms | 2 d / 11 ms | 4 d / 3375 ms | 3 d / 15 ms | FEWER-DIGITS(0<4) |
| bd10.3 machine | 3 | — | 4 d / 404 ms | 4 d / 175 ms | 3 d / 6385 ms | 3 d / 141 ms | FEWER-DIGITS(0<4) |
| bd10.3 wp50 | 2 | — | 2 d / 38 ms | 3 d / 10 ms | 2 d / 3092 ms | 5 d / 14 ms | FEWER-DIGITS(0<5) |
| bd10.4 machine | 3 | — | 6 d / 405 ms | 4 d / 170 ms | 4 d / 6362 ms | 3 d / 107 ms | FEWER-DIGITS(0<6) |
| bd10.4 wp50 | 2 | — | 3 d / 36 ms | 3 d / 7.1 ms | 3 d / 590 ms | 3 d / 8.1 ms | FEWER-DIGITS(0<3) |
| bd10.5 machine | 3 | — | 3 d / 406 ms | 4 d / 174 ms | 4 d / 6387 ms | 3 d / 114 ms | FEWER-DIGITS(0<4) |
| bd10.5 wp50 | 2 | — | 2 d / 32 ms | 3 d / 9.6 ms | 2 d / 103 ms | 3 d / 13 ms | FEWER-DIGITS(0<3) |
| bd10.6 machine | 3 | — | 4 d / 406 ms | 3 d / 174 ms | 4 d / 6391 ms | 3 d / 106 ms | FEWER-DIGITS(0<4) |
| bd10.6 wp50 | 2 | — | 2 d / 47 ms | 2 d / 14 ms | 3 d / 2914 ms | 3 d / 6.9 ms | FEWER-DIGITS(0<3) |
| bd10.7 machine | 3 | — | 5 d / 405 ms | 3 d / 170 ms | 5 d / 6391 ms | 4 d / 107 ms | FEWER-DIGITS(0<5) |
| bd10.7 wp50 | 2 | — | 2 d / 33 ms | 2 d / 7.7 ms | 3 d / 1077 ms | 3 d / 3.5 ms | FEWER-DIGITS(0<3) |
| bd10.8 machine | 3 | — | 4 d / 405 ms | 3 d / 169 ms | 4 d / 6402 ms | 3 d / 106 ms | FEWER-DIGITS(0<4) |
| bd10.8 wp50 | 2 | — | 2 d / 44 ms | 2 d / 14 ms | 3 d / 9321 ms | 3 d / 8.2 ms | FEWER-DIGITS(0<3) |
| bd10.9 machine | 2 | — | 6 d / 406 ms | 3 d / 168 ms | 4 d / 6382 ms | 3 d / 146 ms | FEWER-DIGITS(0<6) |
| bd10.9 wp50 | 1 | — | 2 d / 37 ms | 2 d / 9.7 ms | 3 d / 2477 ms | 3 d / 4.7 ms | FEWER-DIGITS(0<3) |
| bd10.10 machine | 2 | — | 3 d / 408 ms | 4 d / 173 ms | 3 d / 6411 ms | 3 d / 149 ms | FEWER-DIGITS(0<4) |
| bd10.10 wp50 | 1 | — | 2 d / 42 ms | 2 d / 8.9 ms | 2 d / 1769 ms | 3 d / 4.3 ms | FEWER-DIGITS(0<3) |

### 246-siam100b-d11-vector-shooting

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd11.1 machine | 10 | — | 14 d / 49 ms | 14 d / 15 ms | — | — | FEWER-DIGITS(0<14) |
| bd11.1 wp50 | 48 | — | ✗0 d / 7466 ms | 50 d / 5385 ms | — | — | FEWER-DIGITS(0<50) |
| bd11.2 machine | 10 | — | 13 d / 10 ms | 13 d / 12 ms | — | — | FEWER-DIGITS(0<13) |
| bd11.2 wp50 | 48 | — | ~2 d / 3773 ms | 50 d / 3985 ms | — | — | FEWER-DIGITS(0<50) |
| bd11.3 machine | 10 | — | 12 d / 26 ms | 12 d / 22 ms | — | — | FEWER-DIGITS(0<12) |
| bd11.3 wp50 | 48 | — | — | 50 d / 28289 ms | — | — | FEWER-DIGITS(0<50) |
| bd11.4 machine | 10 | — | 12 d / 11 ms | 13 d / 19 ms | — | — | FEWER-DIGITS(0<13) |
| bd11.4 wp50 | 48 | — | — | 50 d / 61012 ms | — | — | FEWER-DIGITS(0<50) |
| bd11.5 machine | 10 | — | 11 d / 6.8 ms | 12 d / 25 ms | — | — | FEWER-DIGITS(0<12) |
| bd11.5 wp50 | 48 | — | ✗0 d / 3939 ms | 50 d / 2090 ms | — | — | FEWER-DIGITS(0<50) |
| bd11.6 machine | 10 | — | 14 d / 11 ms | 13 d / 16 ms | — | — | FEWER-DIGITS(0<14) |
| bd11.6 wp50 | 48 | — | ✗0 d / 8369 ms | 50 d / 84827 ms | — | — | FEWER-DIGITS(0<50) |
| bd11.7 machine | 10 | — | 12 d / 6.9 ms | 13 d / 7.7 ms | — | — | FEWER-DIGITS(0<13) |
| bd11.7 wp50 | 48 | — | 50 d / 42034 ms | 50 d / 6264 ms | — | — | FEWER-DIGITS(0<50) |
| bd11.8 machine | 10 | — | 12 d / 36 ms | 14 d / 85 ms | — | — | FEWER-DIGITS(0<14) |
| bd11.8 wp50 | 48 | — | — | 50 d / 98568 ms | — | — | FEWER-DIGITS(0<50) |
| bd11.9 machine | 10 | — | 12 d / 16 ms | 14 d / 8.5 ms | — | — | FEWER-DIGITS(0<14) |
| bd11.9 wp50 | 48 | — | ~1 d / 4773 ms | 50 d / 5499 ms | — | — | FEWER-DIGITS(0<50) |
| bd11.10 machine | 10 | — | 11 d / 122 ms | 13 d / 122 ms | — | — | FEWER-DIGITS(0<13) |
| bd11.10 wp50 | 48 | — | ~15 d / 27217 ms | — | — | — | FEWER-DIGITS(0<15) |

### 247-siam100b-d12-event-time-root

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd12.1 machine | 10 | — | 10 d / 17 ms | 13 d / 16 ms | 12 d / 232 ms | — | FEWER-DIGITS(0<13) |
| bd12.1 wp50 | 48 | — | ~46 d / 111961 ms | 50 d / 34358 ms | — | — | FEWER-DIGITS(0<50) |
| bd12.2 machine | 10 | — | 10 d / 13 ms | 11 d / 34 ms | 11 d / 131 ms | — | FEWER-DIGITS(0<11) |
| bd12.2 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd12.3 machine | 10 | — | 11 d / 9.2 ms | 12 d / 22 ms | 11 d / 133 ms | — | FEWER-DIGITS(0<12) |
| bd12.3 wp50 | 48 | — | — | 50 d / 7409 ms | — | — | FEWER-DIGITS(0<50) |
| bd12.4 machine | 10 | — | 10 d / 13 ms | 11 d / 19 ms | 11 d / 81 ms | — | FEWER-DIGITS(0<11) |
| bd12.4 wp50 | 48 | — | — | 50 d / 6095 ms | — | — | FEWER-DIGITS(0<50) |
| bd12.5 machine | 10 | — | 10 d / 17 ms | 11 d / 32 ms | 10 d / 249 ms | — | FEWER-DIGITS(0<11) |
| bd12.5 wp50 | 48 | — | — | 50 d / 29546 ms | — | — | FEWER-DIGITS(0<50) |
| bd12.6 machine | 10 | — | 10 d / 12 ms | 11 d / 13 ms | 13 d / 74 ms | — | FEWER-DIGITS(0<13) |
| bd12.6 wp50 | 48 | — | — | 50 d / 2191 ms | — | — | FEWER-DIGITS(0<50) |
| bd12.7 machine | 10 | — | 10 d / 10 ms | 12 d / 27 ms | 12 d / 27 ms | — | FEWER-DIGITS(0<12) |
| bd12.7 wp50 | 48 | — | — | 50 d / 2580 ms | — | — | FEWER-DIGITS(0<50) |
| bd12.8 machine | 10 | — | 10 d / 17 ms | 12 d / 18 ms | 11 d / 354 ms | — | FEWER-DIGITS(0<12) |
| bd12.8 wp50 | 48 | — | ~46 d / 53112 ms | 50 d / 16596 ms | — | — | FEWER-DIGITS(0<50) |
| bd12.9 machine | 10 | — | 11 d / 9.5 ms | 12 d / 18 ms | 13 d / 173 ms | — | FEWER-DIGITS(0<13) |
| bd12.9 wp50 | 48 | — | — | 50 d / 7275 ms | — | — | FEWER-DIGITS(0<50) |
| bd12.10 machine | 10 | — | 11 d / 69 ms | ~4 d / 573 ms | 13 d / 165 ms | — | FEWER-DIGITS(0<13) |
| bd12.10 wp50 | 48 | — | ~46 d / 43057 ms | 50 d / 48308 ms | — | — | FEWER-DIGITS(0<50) |

### 248-siam100b-d13-initial-value-threshold

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd13.1 machine | 10 | — | 11 d / 18 ms | 12 d / 41 ms | 15 d / 346 ms | — | FEWER-DIGITS(0<15) |
| bd13.1 wp50 | 48 | — | — | 50 d / 7504 ms | — | — | FEWER-DIGITS(0<50) |
| bd13.2 machine | 10 | — | 12 d / 3.7 ms | 11 d / 24 ms | 12 d / 416 ms | — | FEWER-DIGITS(0<12) |
| bd13.2 wp50 | 48 | — | — | 50 d / 5680 ms | — | — | FEWER-DIGITS(0<50) |
| bd13.3 machine | 10 | — | 11 d / 4.9 ms | 11 d / 42 ms | 12 d / 71 ms | — | FEWER-DIGITS(0<12) |
| bd13.3 wp50 | 48 | — | — | 50 d / 6121 ms | — | — | FEWER-DIGITS(0<50) |
| bd13.4 machine | 10 | — | 10 d / 112 ms | 11 d / 211 ms | — | — | FEWER-DIGITS(0<11) |
| bd13.4 wp50 | 48 | — | — | 50 d / 79239 ms | — | — | FEWER-DIGITS(0<50) |
| bd13.5 machine | 10 | — | 11 d / 36 ms | 13 d / 201 ms | 13 d / 549 ms | — | FEWER-DIGITS(0<13) |
| bd13.5 wp50 | 48 | — | — | 50 d / 108278 ms | — | — | FEWER-DIGITS(0<50) |
| bd13.6 machine | 10 | — | 11 d / 47 ms | 12 d / 147 ms | 11 d / 107 ms | — | FEWER-DIGITS(0<12) |
| bd13.6 wp50 | 48 | — | ~46 d / 35038 ms | 50 d / 12643 ms | — | — | FEWER-DIGITS(0<50) |
| bd13.7 machine | 10 | — | 10 d / 12 ms | 11 d / 25 ms | 14 d / 292 ms | — | FEWER-DIGITS(0<14) |
| bd13.7 wp50 | 48 | — | — | 50 d / 7582 ms | — | — | FEWER-DIGITS(0<50) |
| bd13.8 machine | 10 | — | 11 d / 8.3 ms | 12 d / 22 ms | 13 d / 290 ms | — | FEWER-DIGITS(0<13) |
| bd13.8 wp50 | 48 | — | — | 50 d / 14995 ms | — | — | FEWER-DIGITS(0<50) |
| bd13.9 machine | 10 | — | 10 d / 48 ms | 14 d / 302 ms | 11 d / 837 ms | — | FEWER-DIGITS(0<14) |
| bd13.9 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd13.10 machine | 10 | — | 11 d / 9.4 ms | 13 d / 62 ms | 12 d / 280 ms | — | FEWER-DIGITS(0<13) |
| bd13.10 wp50 | 48 | — | ~47 d / 45046 ms | 50 d / 16194 ms | — | — | FEWER-DIGITS(0<50) |

### 249-siam100b-d14-fem-parameter

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd14.1 machine | 10 | — | ~9 d / 11938 ms | — | — | — | FEWER-DIGITS(0<9) |
| bd14.1 wp50 | 48 | — | ~9 d / 24025 ms | — | — | — | FEWER-DIGITS(0<9) |
| bd14.2 machine | 10 | — | ~8 d / 7683 ms | — | — | — | FEWER-DIGITS(0<8) |
| bd14.2 wp50 | 48 | — | ~10 d / 115058 ms | — | — | — | FEWER-DIGITS(0<10) |
| bd14.3 machine | 10 | — | ~7 d / 5272 ms | — | — | — | FEWER-DIGITS(0<7) |
| bd14.3 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd14.4 machine | 10 | — | ~8 d / 11002 ms | — | — | — | FEWER-DIGITS(0<8) |
| bd14.4 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd14.5 machine | 10 | — | ~9 d / 40218 ms | — | — | — | FEWER-DIGITS(0<9) |
| bd14.5 wp50 | 48 | — | ~9 d / 57792 ms | — | — | — | FEWER-DIGITS(0<9) |
| bd14.6 machine | 10 | — | ~8 d / 9635 ms | — | — | — | FEWER-DIGITS(0<8) |
| bd14.6 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd14.7 machine | 10 | — | ~8 d / 5736 ms | — | — | — | FEWER-DIGITS(0<8) |
| bd14.7 wp50 | 48 | — | ~8 d / 57314 ms | — | — | — | FEWER-DIGITS(0<8) |
| bd14.8 machine | 10 | — | ~1 d / 5439 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd14.8 wp50 | 48 | — | ~1 d / 11602 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd14.9 machine | 10 | — | ~7 d / 4028 ms | — | — | — | FEWER-DIGITS(0<7) |
| bd14.9 wp50 | 48 | — | ~8 d / 116877 ms | — | — | — | FEWER-DIGITS(0<8) |
| bd14.10 machine | 10 | — | ~1 d / 21303 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd14.10 wp50 | 48 | — | ~1 d / 30014 ms | — | — | — | FEWER-DIGITS(0<1) |

### 250-siam100b-d15-blowup-times

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd15.1 machine | 10 | — | ~2 d / 46 ms | ~7 d / 581 ms | — | — | FEWER-DIGITS(0<7) |
| bd15.1 wp50 | 48 | — | ~16 d / 14888 ms | — | — | — | FEWER-DIGITS(0<16) |
| bd15.2 machine | 10 | — | ~2 d / 62 ms | ~2 d / 349 ms | — | — | FEWER-DIGITS(0<2) |
| bd15.2 wp50 | 48 | — | ~2 d / 6862 ms | — | — | — | FEWER-DIGITS(0<2) |
| bd15.3 machine | 10 | — | ~3 d / 157 ms | ~4 d / 1461 ms | — | — | FEWER-DIGITS(0<4) |
| bd15.3 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd15.4 machine | 10 | — | ~2 d / 17 ms | ~3 d / 453 ms | — | — | FEWER-DIGITS(0<3) |
| bd15.4 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd15.5 machine | 10 | — | ~2 d / 2318 ms | ~7 d / 152 ms | — | — | FEWER-DIGITS(0<7) |
| bd15.5 wp50 | 48 | — | ~15 d / 66145 ms | — | — | — | FEWER-DIGITS(0<15) |
| bd15.6 machine | 10 | — | ~6 d / 2.8 ms | 13 d / 55 ms | — | — | FEWER-DIGITS(0<13) |
| bd15.6 wp50 | 48 | — | ✗0 d / 50 ms | — | — | — | OK(nobody) |
| bd15.7 machine | 10 | — | ~6 d / 2.6 ms | 13 d / 88 ms | — | — | FEWER-DIGITS(0<13) |
| bd15.7 wp50 | 48 | — | ~2 d / 613 ms | — | — | — | FEWER-DIGITS(0<2) |
| bd15.8 machine | 10 | — | ~2 d / 36 ms | ~7 d / 626 ms | — | — | FEWER-DIGITS(0<7) |
| bd15.8 wp50 | 48 | — | ~16 d / 51405 ms | — | — | — | FEWER-DIGITS(0<16) |
| bd15.9 machine | 10 | — | ~3 d / 4351 ms | ~4 d / 380 ms | — | — | FEWER-DIGITS(0<4) |
| bd15.9 wp50 | 48 | — | ✗0 d / 178 ms | — | — | — | OK(nobody) |
| bd15.10 machine | 10 | — | ~1 d / 245 ms | ~1 d / 2158 ms | — | — | FEWER-DIGITS(0<1) |
| bd15.10 wp50 | 48 | — | ~1 d / 104669 ms | — | — | — | FEWER-DIGITS(0<1) |

### 251-siam100b-d16-heat-functionals

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd16.1 machine | 10 | — | ✗0 d / 89 ms | ~6 d / 226 ms | — | — | FEWER-DIGITS(0<6) |
| bd16.1 wp50 | 48 | — | ✗0 d / 6628 ms | — | — | — | OK(nobody) |
| bd16.2 machine | 10 | — | ~4 d / 2774 ms | ~5 d / 107 ms | — | — | FEWER-DIGITS(0<5) |
| bd16.2 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd16.3 machine | 10 | — | ~4 d / 5649 ms | ~4 d / 186 ms | — | — | FEWER-DIGITS(0<4) |
| bd16.3 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd16.4 machine | 10 | — | ✗0 d / 510 ms | ~3 d / 133 ms | — | — | FEWER-DIGITS(0<3) |
| bd16.4 wp50 | 48 | — | ✗0 d / 41232 ms | — | — | — | OK(nobody) |
| bd16.5 machine | 10 | — | ~3 d / 1170 ms | ~6 d / 208 ms | — | — | FEWER-DIGITS(0<6) |
| bd16.5 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd16.6 machine | 10 | — | ~3 d / 145 ms | ~5 d / 103 ms | — | — | FEWER-DIGITS(0<5) |
| bd16.6 wp50 | 48 | — | ~7 d / 5915 ms | — | — | — | FEWER-DIGITS(0<7) |
| bd16.7 machine | 10 | — | ~2 d / 5.1 ms | ~5 d / 99 ms | — | — | FEWER-DIGITS(0<5) |
| bd16.7 wp50 | 48 | — | ~6 d / 6589 ms | — | — | — | FEWER-DIGITS(0<6) |
| bd16.8 machine | 10 | — | ✗0 d / 66 ms | ~6 d / 223 ms | — | — | FEWER-DIGITS(0<6) |
| bd16.8 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd16.9 machine | 10 | — | ✗0 d / 45 ms | ~4 d / 4061 ms | — | — | FEWER-DIGITS(0<4) |
| bd16.9 wp50 | 48 | — | ✗0 d / 39205 ms | — | — | — | OK(nobody) |
| bd16.10 machine | 10 | — | ~3 d / 3783 ms | ~2 d / 113 ms | — | — | FEWER-DIGITS(0<3) |
| bd16.10 wp50 | 48 | — | — | — | — | — | OK(nobody) |

### 252-siam100b-d17-nodal-lines

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd17.1 machine | 10 | — | 14 d / 2503 ms | — | — | — | FEWER-DIGITS(0<14) |
| bd17.1 wp50 | 48 | — | ~14 d / 504 ms | — | — | — | FEWER-DIGITS(0<14) |
| bd17.2 machine | 10 | — | ~1 d / 648 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.2 wp50 | 48 | — | ~1 d / 661 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.3 machine | 10 | — | ~4 d / 549 ms | — | — | — | FEWER-DIGITS(0<4) |
| bd17.3 wp50 | 48 | — | ~4 d / 472 ms | — | — | — | FEWER-DIGITS(0<4) |
| bd17.4 machine | 10 | — | ~1 d / 944 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.4 wp50 | 48 | — | ~1 d / 913 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.5 machine | 10 | — | ~1 d / 652 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.5 wp50 | 48 | — | ~1 d / 642 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.6 machine | 10 | — | ~1 d / 576 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.6 wp50 | 48 | — | ~1 d / 577 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.7 machine | 10 | — | ~1 d / 701 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.7 wp50 | 48 | — | ~1 d / 696 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.8 machine | 10 | — | ~1 d / 368 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.8 wp50 | 48 | — | ~1 d / 368 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd17.9 machine | 10 | — | ~2 d / 961 ms | — | — | — | FEWER-DIGITS(0<2) |
| bd17.9 wp50 | 48 | — | ~2 d / 918 ms | — | — | — | FEWER-DIGITS(0<2) |
| bd17.10 machine | 10 | — | ~3 d / 581 ms | — | — | — | FEWER-DIGITS(0<3) |
| bd17.10 wp50 | 48 | — | ~3 d / 585 ms | — | — | — | FEWER-DIGITS(0<3) |

### 253-siam100b-d18-dirichlet-eigenvalue

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd18.1 machine | 10 | — | ~6 d / 137 ms | — | — | — | FEWER-DIGITS(0<6) |
| bd18.1 wp50 | 48 | — | ~6 d / 154 ms | — | — | — | FEWER-DIGITS(0<6) |
| bd18.2 machine | 10 | — | ~4 d / 480 ms | — | — | — | FEWER-DIGITS(0<4) |
| bd18.2 wp50 | 48 | — | ~4 d / 45 ms | — | — | — | FEWER-DIGITS(0<4) |
| bd18.3 machine | 10 | — | ~4 d / 45 ms | — | — | — | FEWER-DIGITS(0<4) |
| bd18.3 wp50 | 48 | — | ~4 d / 45 ms | — | — | — | FEWER-DIGITS(0<4) |
| bd18.4 machine | 10 | — | ~5 d / 95 ms | — | — | — | FEWER-DIGITS(0<5) |
| bd18.4 wp50 | 48 | — | ~5 d / 90 ms | — | — | — | FEWER-DIGITS(0<5) |
| bd18.5 machine | 10 | — | ~2 d / 95 ms | — | — | — | FEWER-DIGITS(0<2) |
| bd18.5 wp50 | 48 | — | ~2 d / 86 ms | — | — | — | FEWER-DIGITS(0<2) |
| bd18.6 machine | 10 | — | ~3 d / 84 ms | — | — | — | FEWER-DIGITS(0<3) |
| bd18.6 wp50 | 48 | — | ~3 d / 83 ms | — | — | — | FEWER-DIGITS(0<3) |
| bd18.7 machine | 10 | — | ~6 d / 11324 ms | — | — | — | FEWER-DIGITS(0<6) |
| bd18.7 wp50 | 48 | — | ~6 d / 11966 ms | — | — | — | FEWER-DIGITS(0<6) |
| bd18.8 machine | 10 | — | ~1 d / 2074 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd18.8 wp50 | 48 | — | ~1 d / 1686 ms | — | — | — | FEWER-DIGITS(0<1) |
| bd18.9 machine | 10 | — | ~2 d / 3472 ms | — | — | — | FEWER-DIGITS(0<2) |
| bd18.9 wp50 | 48 | — | ~2 d / 2654 ms | — | — | — | FEWER-DIGITS(0<2) |
| bd18.10 machine | 10 | — | ~2 d / 110 ms | — | — | — | FEWER-DIGITS(0<2) |
| bd18.10 wp50 | 48 | — | ~2 d / 113 ms | — | — | — | FEWER-DIGITS(0<2) |

### 254-siam100b-d19-design

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd19.1 machine | 10 | — | ~5 d / 522 ms | 13 d / 128 ms | 14 d / 218 ms | — | FEWER-DIGITS(0<14) |
| bd19.1 wp50 | 48 | — | ~9 d / 692 ms | — | — | — | FEWER-DIGITS(0<9) |
| bd19.2 machine | 10 | — | ~3 d / 736 ms | ~5 d / 342 ms | 15 d / 487 ms | — | FEWER-DIGITS(0<15) |
| bd19.2 wp50 | 48 | — | ~10 d / 1937 ms | — | — | — | FEWER-DIGITS(0<10) |
| bd19.3 machine | 10 | — | ~4 d / 2781 ms | 13 d / 134 ms | 14 d / 2197 ms | — | FEWER-DIGITS(0<14) |
| bd19.3 wp50 | 48 | — | ~9 d / 8047 ms | — | — | — | FEWER-DIGITS(0<9) |
| bd19.4 machine | 10 | — | ~2 d / 25933 ms | ~4 d / 5787 ms | 12 d / 26954 ms | — | FEWER-DIGITS(0<12) |
| bd19.4 wp50 | 48 | — | ~9 d / 72020 ms | — | — | — | FEWER-DIGITS(0<9) |
| bd19.5 machine | 10 | — | ~3 d / 39504 ms | ~7 d / 5550 ms | — | — | FEWER-DIGITS(0<7) |
| bd19.5 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd19.6 machine | 10 | — | ~5 d / 71752 ms | ~6 d / 11281 ms | — | — | FEWER-DIGITS(0<6) |
| bd19.6 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| bd19.7 machine | 10 | — | ~5 d / 5437 ms | 14 d / 708 ms | 15 d / 1529 ms | — | FEWER-DIGITS(0<15) |
| bd19.7 wp50 | 48 | — | ~8 d / 13807 ms | — | — | — | FEWER-DIGITS(0<8) |
| bd19.8 machine | 10 | — | ~2 d / 2625 ms | 10 d / 356 ms | 13 d / 2769 ms | — | FEWER-DIGITS(0<13) |
| bd19.8 wp50 | 48 | — | ~11 d / 7751 ms | — | — | — | FEWER-DIGITS(0<11) |
| bd19.9 machine | 10 | — | ~1 d / 10491 ms | ~3 d / 3827 ms | 13 d / 13727 ms | — | FEWER-DIGITS(0<13) |
| bd19.9 wp50 | 48 | — | ~8 d / 34863 ms | — | — | — | FEWER-DIGITS(0<8) |
| bd19.10 machine | 10 | — | — | — | — | — | OK(nobody) |
| bd19.10 wp50 | 48 | — | — | — | — | — | OK(nobody) |

### 255-siam100b-d20-functional-roots

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| bd20.1 machine | 10 | — | 15 d / 0.56 ms | 15 d / 0.099 ms | 15 d / 4.4 ms | 15 d / 163 ms | FEWER-DIGITS(0<15) |
| bd20.1 wp50 | 48 | — | 50 d / 2.2 ms | 50 d / 4.4 ms | 50 d / 4.8 ms | 50 d / 13 ms | FEWER-DIGITS(0<50) |
| bd20.2 machine | 10 | — | 15 d / 2.5 ms | 15 d / 0.22 ms | 15 d / 24 ms | 15 d / 36 ms | FEWER-DIGITS(0<15) |
| bd20.2 wp50 | 48 | — | 50 d / 5 ms | 50 d / 21 ms | 50 d / 199 ms | 50 d / 32 ms | FEWER-DIGITS(0<50) |
| bd20.3 machine | 10 | — | 15 d / 7.8 ms | 15 d / 0.46 ms | 15 d / 32 ms | 15 d / 24 ms | FEWER-DIGITS(0<15) |
| bd20.3 wp50 | 48 | — | 50 d / 12 ms | 50 d / 61 ms | 50 d / 649 ms | 50 d / 50 ms | FEWER-DIGITS(0<50) |
| bd20.4 machine | 10 | — | ~7 d / 1.3 ms | ~7 d / 0.16 ms | ~7 d / 8 ms | 13 d / 4.9 ms | FEWER-DIGITS(0<13) |
| bd20.4 wp50 | 48 | — | 50 d / 3.3 ms | 50 d / 9.9 ms | 50 d / 69 ms | 50 d / 8.3 ms | FEWER-DIGITS(0<50) |
| bd20.5 machine | 10 | — | 13 d / 5.2 ms | 13 d / 0.35 ms | 15 d / 33 ms | 15 d / 246 ms | FEWER-DIGITS(0<15) |
| bd20.5 wp50 | 48 | — | 50 d / 14 ms | 50 d / 56 ms | 50 d / 26 ms | 50 d / 249 ms | FEWER-DIGITS(0<50) |
| bd20.6 machine | 10 | — | 13 d / 2.2 ms | 13 d / 0.17 ms | 15 d / 10 ms | 15 d / 100 ms | FEWER-DIGITS(0<15) |
| bd20.6 wp50 | 48 | — | 50 d / 5.1 ms | 50 d / 31 ms | 50 d / 14 ms | 50 d / 100 ms | FEWER-DIGITS(0<50) |
| bd20.7 machine | 10 | — | 15 d / 1.3 ms | 15 d / 0.12 ms | — | — | FEWER-DIGITS(0<15) |
| bd20.7 wp50 | 48 | — | 50 d / 24 ms | 50 d / 2907 ms | — | — | FEWER-DIGITS(0<50) |
| bd20.8 machine | 10 | — | 15 d / 0.6 ms | 15 d / 2 ms | — | — | FEWER-DIGITS(0<15) |
| bd20.8 wp50 | 48 | — | 50 d / 9.6 ms | 50 d / 142 ms | — | — | FEWER-DIGITS(0<50) |
| bd20.9 machine | 10 | — | 15 d / 2.4 ms | 15 d / 0.17 ms | — | — | FEWER-DIGITS(0<15) |
| bd20.9 wp50 | 48 | — | 50 d / 34 ms | 50 d / 5106 ms | — | — | FEWER-DIGITS(0<50) |
| bd20.10 machine | 7 | — | 8 d / 18 ms | 8 d / 0.85 ms | 15 d / 131 ms | 15 d / 1251 ms | FEWER-DIGITS(0<15) |
| bd20.10 wp50 | 48 | — | 50 d / 52 ms | 50 d / 174 ms | 50 d / 97 ms | 50 d / 1358 ms | FEWER-DIGITS(0<50) |
