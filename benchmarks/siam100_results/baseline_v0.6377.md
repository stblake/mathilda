# SIAM 100-digit scoring

Cell: correct digits / ms (`~` below need, `✗` wrong, — none).

| experiment | cases | Mathilda on target |
|---|---|---|
| 166-siam100-c1-twisted-tail | 2 | 1 |
| 167-siam100-c2-photon-mirrors | 3 | 1 |
| 168-siam100-c3-operator-norm | 4 | 2 |
| 169-siam100-c4-global-minimum | 2 | 0 |
| 170-siam100-c5-gamma-best-cubic | 2 | 1 |
| 171-siam100-c6-biased-flea | 2 | 1 |
| 172-siam100-c7-sparse-inverse | 4 | 3 |
| 173-siam100-c8-heat-plate | 3 | 1 |
| 174-siam100-c9-parametric-maximum | 2 | 0 |
| 175-siam100-c10-brownian-ends | 3 | 1 |
| 176-siam100-d1-kempner-42 | 2 | 2 |
| 177-siam100-d2-digit-length-series | 2 | 2 |
| 178-siam100-d3-thomson-100 | 2 | 2 |
| 179-siam100-d4-riemann-r-zero | 3 | 3 |
| 180-siam100-d5-toeplitz-polynomial | 2 | 2 |
| 181-siam100-d6-quadruple-exponential | 2 | 1 |
| 182-siam100-d7-region-area | 2 | 0 |
| 183-siam100-d8-sobolev-constant | 4 | 4 |
| 184-siam100-d9-triangular-walk | 2 | 1 |
| 185-siam100-d10-random-recurrence | 2 | 0 |
| 186-siam100-d11-six-masses | 2 | 1 |
| 187-siam100-d12-constrained-particle | 2 | 1 |
| 188-siam100-d13-nonlinear-ode-z0 | 2 | 0 |
| 189-siam100-d14-pentagon-poisson | 2 | 1 |
| 190-siam100-d15-blowup-time | 2 | 1 |
| 191-siam100-d16-daubechies-heat | 3 | 3 |
| 192-siam100-d17-lshape-nodal-line | 3 | 3 |
| 193-siam100-d18-koch-eigenvalue | 3 | 3 |
| 194-siam100-d19-optimal-design | 3 | 1 |
| 195-siam100-d20-half-exponential | 4 | 2 |

| system | correct | inexact | wrong | none |
|---|---|---|---|---|
| Mathilda | 25 | 9 | 8 | 34 |
| WL | 21 | 10 | 5 | 40 |
| Py | 33 | 6 | 7 | 30 |
| Maxima | 20 | 2 | 1 | 53 |
| FriCAS | 14 | 0 | 2 | 60 |

### 166-siam100-c1-twisted-tail

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c1.machine | 10 | ~9 d / 29 ms | ~8 d / 25 ms | ✗0 d / 11 ms | ✗0 d / 201 ms | ✗0 d / 446 ms | OK |
| c1.wp50 | 48 | ✗0 d / 16115 ms | ~9 d / 219 ms | ✗0 d / 95 ms | — | ✗0 d / 109454 ms | FEWER-DIGITS(0<9) |

### 167-siam100-c2-photon-mirrors

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c2.machine | 3 | 5 d / 12 ms | 4 d / 42 ms | 4 d / 4 ms | 5 d / 13 ms | 5 d / 140 ms | SLOWER(2.9x) |
| c2.machine-t3 | 10 | 14 d / 1.9 ms | 13 d / 7.3 ms | 13 d / 0.72 ms | 14 d / 2 ms | 14 d / 7.8 ms | SLOWER(2.6x) |
| c2.wp50 | 45 | 50 d / 29 ms | 45 d / 97 ms | 50 d / 50 ms | 49 d / 200 ms | 49 d / 244 ms | OK |

### 168-siam100-c3-operator-norm

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c3.machine | 10 | 10 d / 57063 ms | 10 d / 46260 ms | 10 d / 28910 ms | — | — | SLOWER(2.0x) |
| c3.wp50 | 38 | — | — | — | — | — | OK(nobody) |
| c3.machine-n500 | 10 | 15 d / 447 ms | 15 d / 717 ms | 15 d / 49 ms | 15 d / 3080 ms | — | SLOWER(9.2x) |
| c3.wp50-n100 | 48 | 50 d / 264 ms | 50 d / 578 ms | 50 d / 7163 ms | — | — | OK |

### 169-siam100-c4-global-minimum

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c4.machine | 10 | 14 d / 12 ms | 15 d / 240 ms | 15 d / 147 ms | — | — | FEWER-DIGITS(14<15) |
| c4.wp50 | 48 | ~15 d / 623 ms | ✗0 d / 1763 ms | 50 d / 182 ms | — | — | FEWER-DIGITS(15<50) |

### 170-siam100-c5-gamma-best-cubic

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c5.machine | 10 | ~7 d / 1551 ms | — | 11 d / 19323 ms | — | — | FEWER-DIGITS(7<11) |
| c5.wp50 | 48 | — | — | — | — | — | OK(nobody) |

### 171-siam100-c6-biased-flea

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c6.machine | 10 | ~7 d / 44 ms | — | 12 d / 844 ms | 12 d / 39625 ms | — | FEWER-DIGITS(7<12) |
| c6.wp50 | 48 | — | — | — | — | — | OK(nobody) |

### 172-siam100-c7-sparse-inverse

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c7.machine | 10 | 12 d / 227 ms | — | 15 d / 2679 ms | — | — | FEWER-DIGITS(12<15) |
| c7.wp50 | 48 | — | — | — | — | — | OK(nobody) |
| c7.machine n=2000 | 10 | 15 d / 89 ms | — | 15 d / 1069 ms | — | 14 d / 160771 ms | OK |
| c7.wp50 n=2000 | 48 | 50 d / 41493 ms | — | — | — | — | OK |

### 173-siam100-c8-heat-plate

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c8.machine | 10 | 15 d / 11 ms | 15 d / 50 ms | 15 d / 5.8 ms | 15 d / 13 ms | — | SLOWER(1.8x) |
| c8.wp50 | 48 | 50 d / 2354 ms | 50 d / 3755 ms | 50 d / 10 ms | 50 d / 98 ms | — | SLOWER(229.8x) |
| c8.machine pde | 4 | — | ✗0 d / 760 ms | — | — | — | OK(nobody) |

### 174-siam100-c9-parametric-maximum

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c9.machine | 10 | ~2 d / 15503 ms | ~6 d / 98539 ms | ~3 d / 35126 ms | — | — | FEWER-DIGITS(2<6) |
| c9.wp50 | 48 | — | ~30 d / 9224 ms | — | — | — | FEWER-DIGITS(0<30) |

### 175-siam100-c10-brownian-ends

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c10.machine | 10 | 15 d / 4.2 ms | — | 15 d / 0.55 ms | 15 d / 0.85 ms | 15 d / 25 ms | SLOWER(7.7x) |
| c10.wp50 | 48 | — | — | 50 d / 1.1 ms | 50 d / 4 ms | 50 d / 25 ms | FEWER-DIGITS(0<50) |
| c10.machine pde | 3 | — | — | — | — | — | OK(nobody) |

### 176-siam100-d1-kempner-42

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d1.machine | 10 | ✗0 d / 1879 ms | ✗0 d / 606 ms | ✗0 d / 72 ms | — | — | OK(nobody) |
| d1.wp50 | 48 | — | — | ✗0 d / 1660 ms | — | — | OK(nobody) |

### 177-siam100-d2-digit-length-series

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d2.machine | 10 | ~1 d / 2.4 ms | ~1 d / 27 ms | ~1 d / 59 ms | — | — | OK |
| d2.wp50 | 48 | ~1 d / 12 ms | ~1 d / 31 ms | ~1 d / 1614 ms | — | — | OK |

### 178-siam100-d3-thomson-100

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d3.machine | 10 | 12 d / 15980 ms | — | ~9 d / 45152 ms | — | — | OK |
| d3.wp50 | 48 | — | — | — | — | — | OK(nobody) |

### 179-siam100-d4-riemann-r-zero

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d4.machine | 10 | ✗0 d / 5.7 ms | — | ✗0 d / 64 ms | — | — | OK(nobody) |
| d4.wp50 | 48 | ✗0 d / 25709 ms | — | ✗0 d / 112 ms | — | — | OK(nobody) |
| d4.wp50-builtin | 48 | — | — | — | — | — | OK(nobody) |

### 180-siam100-d5-toeplitz-polynomial

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d5.machine | 10 | ✗0 d / 4241 ms | ✗0 d / 7847 ms | ✗0 d / 866 ms | — | — | OK(nobody) |
| d5.wp50 | 48 | — | ✗0 d / 581693 ms | — | — | — | OK(nobody) |

### 181-siam100-d6-quadruple-exponential

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d6.machine | 10 | 13 d / 0.13 ms | 12 d / 8.5 ms | 15 d / 0.54 ms | 15 d / 11 ms | 15 d / 105 ms | FEWER-DIGITS(13<15) |
| d6.wp50 | 48 | 50 d / 151 ms | 50 d / 388 ms | 50 d / 265 ms | ~8 d / 688 ms | 50 d / 152660 ms | OK |

### 182-siam100-d7-region-area

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d7.machine | 10 | 11 d / 33 ms | 14 d / 5034 ms | 14 d / 9.7 ms | 15 d / 212 ms | — | FEWER-DIGITS(11<15) SLOWER(3.4x) |
| d7.wp50 | 48 | 50 d / 878 ms | 50 d / 5812 ms | 50 d / 180 ms | ~4 d / 2032 ms | — | SLOWER(4.9x) |

### 183-siam100-d8-sobolev-constant

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d8.machine | 3 | — | — | — | — | — | OK(nobody) |
| d8.wp50 | 3 | — | — | — | — | — | OK(nobody) |
| d8.machine c3R1 | 10 | 15 d / 0.079 ms | 12 d / 3.6 ms | 15 d / 0.14 ms | 15 d / 7 ms | — | OK |
| d8.machine c2R2 | 10 | 15 d / 119 ms | ~9 d / 39 ms | 15 d / 84 ms | 15 d / 1692 ms | — | OK |

### 184-siam100-d9-triangular-walk

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d9.machine | 10 | 14 d / 60 ms | 14 d / 26 ms | 14 d / 15 ms | — | 15 d / 50158 ms | FEWER-DIGITS(14<15) SLOWER(4.1x) |
| d9.wp50 | 48 | 50 d / 1323 ms | 50 d / 1831 ms | 50 d / 13765 ms | 50 d / 3509 ms | 50 d / 2103 ms | OK |

### 185-siam100-d10-random-recurrence

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d10.machine | 2 | 4 d / 1677 ms | 4 d / 187 ms | 4 d / 242 ms | 4 d / 11807 ms | 3 d / 468 ms | SLOWER(9.0x) |
| d10.wp50 | 1 | 2 d / 47 ms | 2 d / 46 ms | 3 d / 11 ms | 2 d / 1812 ms | 3 d / 35 ms | FEWER-DIGITS(2<3) SLOWER(4.5x) |

### 186-siam100-d11-six-masses

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d11.machine | 10 | — | 13 d / 283 ms | 14 d / 1024 ms | — | — | FEWER-DIGITS(0<14) |
| d11.wp50 | 48 | — | — | — | — | — | OK(nobody) |

### 187-siam100-d12-constrained-particle

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d12.machine | 10 | — | 10 d / 276 ms | 12 d / 1899 ms | 11 d / 2409 ms | — | FEWER-DIGITS(0<12) |
| d12.wp50 | 48 | — | — | — | — | — | OK(nobody) |

### 188-siam100-d13-nonlinear-ode-z0

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d13.machine | 10 | — | 10 d / 227 ms | 11 d / 400 ms | — | — | FEWER-DIGITS(0<11) |
| d13.wp50 | 48 | — | — | 50 d / 178353 ms | — | — | FEWER-DIGITS(0<50) |

### 189-siam100-d14-pentagon-poisson

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d14.machine | 10 | — | ~3 d / 16308 ms | — | — | — | FEWER-DIGITS(0<3) |
| d14.wp50 | 10 | — | — | — | — | — | OK(nobody) |

### 190-siam100-d15-blowup-time

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d15.machine | 6 | ~3 d / 4809 ms | ~2 d / 249 ms | ~4 d / 4811 ms | — | — | FEWER-DIGITS(3<4) |
| d15.wp50 | 6 | ~3 d / 4402 ms | ~3 d / 16628 ms | — | — | — | OK |

### 191-siam100-d16-daubechies-heat

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d16.machine | 10 | — | — | — | — | — | OK(nobody) |
| d16.machine-reduced | 10 | — | — | — | — | — | OK(nobody) |
| d16.wp50 | 48 | — | — | — | — | — | OK(nobody) |

### 192-siam100-d17-lshape-nodal-line

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d17.machine | 10 | — | — | — | — | — | OK(nobody) |
| d17.machine-reduced | 10 | — | — | — | — | — | OK(nobody) |
| d17.wp50 | 48 | — | — | — | — | — | OK(nobody) |

### 193-siam100-d18-koch-eigenvalue

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d18.machine | 5 | — | — | — | — | — | OK(nobody) |
| d18.machine level=0 | 10 | — | — | — | — | — | OK(nobody) |
| d18.wp50 | 5 | — | — | — | — | — | OK(nobody) |

### 194-siam100-d19-optimal-design

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d19.machine | 10 | ✗0 d / 1689 ms | — | ~3 d / 28742 ms | — | — | FEWER-DIGITS(0<3) |
| d19.machine-reduced | 10 | ✗0 d / 77 ms | — | 11 d / 209 ms | 14 d / 5603 ms | — | FEWER-DIGITS(0<14) |
| d19.wp50 | 48 | ✗0 d / 62294 ms | — | — | — | — | OK(nobody) |

### 195-siam100-d20-half-exponential

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d20.machine | 5 | — | — | — | — | — | OK(nobody) |
| d20.wp50 | 5 | — | — | — | — | — | OK(nobody) |
| d20.machine robbins n=20 | 10 | 14 d / 2 ms | — | 15 d / 0.5 ms | 15 d / 19 ms | 15 d / 440 ms | FEWER-DIGITS(14<15) SLOWER(3.9x) |
| d20.wp50 robbins n=60 | 48 | — | — | 48 d / 852 ms | 50 d / 745 ms | 50 d / 10368 ms | FEWER-DIGITS(0<50) |
