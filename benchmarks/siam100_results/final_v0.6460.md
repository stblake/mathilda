# SIAM 100-digit scoring

Cell: correct digits / ms (`~` below need, `✗` wrong, — none).

| experiment | cases | Mathilda on target |
|---|---|---|
| 166-siam100-c1-twisted-tail | 2 | 2 |
| 167-siam100-c2-photon-mirrors | 3 | 2 |
| 168-siam100-c3-operator-norm | 4 | 4 |
| 169-siam100-c4-global-minimum | 2 | 2 |
| 170-siam100-c5-gamma-best-cubic | 2 | 2 |
| 171-siam100-c6-biased-flea | 2 | 1 |
| 172-siam100-c7-sparse-inverse | 4 | 4 |
| 173-siam100-c8-heat-plate | 3 | 3 |
| 174-siam100-c9-parametric-maximum | 2 | 2 |
| 175-siam100-c10-brownian-ends | 3 | 3 |
| 176-siam100-d1-kempner-42 | 2 | 2 |
| 177-siam100-d2-digit-length-series | 2 | 2 |
| 178-siam100-d3-thomson-100 | 2 | 2 |
| 179-siam100-d4-riemann-r-zero | 3 | 3 |
| 180-siam100-d5-toeplitz-polynomial | 2 | 2 |
| 181-siam100-d6-quadruple-exponential | 2 | 2 |
| 182-siam100-d7-region-area | 2 | 2 |
| 183-siam100-d8-sobolev-constant | 4 | 4 |
| 184-siam100-d9-triangular-walk | 2 | 1 |
| 185-siam100-d10-random-recurrence | 2 | 1 |
| 186-siam100-d11-six-masses | 2 | 1 |
| 187-siam100-d12-constrained-particle | 2 | 2 |
| 188-siam100-d13-nonlinear-ode-z0 | 2 | 2 |
| 189-siam100-d14-pentagon-poisson | 2 | 2 |
| 190-siam100-d15-blowup-time | 2 | 1 |
| 191-siam100-d16-daubechies-heat | 3 | 3 |
| 192-siam100-d17-lshape-nodal-line | 3 | 3 |
| 193-siam100-d18-koch-eigenvalue | 3 | 3 |
| 194-siam100-d19-optimal-design | 3 | 3 |
| 195-siam100-d20-half-exponential | 4 | 3 |
| 196-siam100h-c1-oscillatory-tails | 20 | 20 |
| 197-siam100h-c2-mirror-billiards | 20 | 13 |
| 198-siam100h-c3-matrix-norms | 20 | 18 |
| 199-siam100h-c4-global-minima | 20 | 20 |
| 200-siam100h-c5-discrete-minimax | 20 | 19 |
| 201-siam100h-c6-lattice-return | 20 | 13 |
| 202-siam100h-c7-sparse-solve | 20 | 20 |
| 203-siam100h-c8-heat-crossing | 20 | 20 |
| 204-siam100h-c9-parametric-extremum | 20 | 19 |
| 205-siam100h-c10-brownian-hitting | 20 | 20 |
| 206-siam100h-d1-kempner-family | 20 | 17 |
| 207-siam100h-d2-digit-length-family | 20 | 16 |
| 208-siam100h-d3-point-energy-family | 20 | 19 |
| 209-siam100h-d4-number-theory-roots-family | 20 | 19 |
| 210-siam100h-d5-matrix-polynomial-family | 20 | 17 |
| 211-siam100h-d6-oscillatory-towers | 20 | 14 |
| 212-siam100h-d7-region-areas | 20 | 19 |
| 213-siam100h-d8-sobolev-family | 20 | 17 |
| 214-siam100h-d9-lattice-walks | 20 | 16 |
| 215-siam100h-d10-random-recurrences | 20 | 6 |
| 216-siam100h-d11-shooting-targets | 20 | 19 |
| 217-siam100h-d12-event-parameter | 20 | 18 |
| 218-siam100h-d13-initial-value-shooting | 20 | 15 |
| 219-siam100h-d14-elliptic-parameter | 20 | 16 |
| 220-siam100h-d15-blowup-times | 20 | 11 |
| 221-siam100h-d16-heat-crossings | 20 | 20 |
| 222-siam100h-d17-nodal-lines | 20 | 20 |
| 223-siam100h-d18-dirichlet-eigenvalues | 20 | 20 |
| 224-siam100h-d19-design-variants | 20 | 15 |
| 225-siam100h-d20-half-iterates | 20 | 10 |

| system | correct | inexact | wrong | none |
|---|---|---|---|---|
| Mathilda | 484 | 127 | 26 | 39 |
| WL | 286 | 185 | 70 | 135 |
| Py | 363 | 64 | 51 | 198 |
| Maxima | 214 | 25 | 6 | 431 |
| FriCAS | 127 | 14 | 14 | 521 |

### 166-siam100-c1-twisted-tail

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c1.machine | 10 | 12 d / 0.2 ms | ~8 d / 25 ms | ✗0 d / 11 ms | ✗0 d / 201 ms | ✗0 d / 446 ms | OK |
| c1.wp50 | 48 | 50 d / 151 ms | ~9 d / 219 ms | ✗0 d / 95 ms | — | ✗0 d / 109454 ms | OK |

### 167-siam100-c2-photon-mirrors

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c2.machine | 3 | 5 d / 6.2 ms | 4 d / 42 ms | 4 d / 4 ms | 5 d / 13 ms | 5 d / 140 ms | SLOWER(1.5x) |
| c2.machine-t3 | 10 | 14 d / 1.1 ms | 13 d / 7.3 ms | 13 d / 0.72 ms | 14 d / 2 ms | 14 d / 7.8 ms | OK |
| c2.wp50 | 45 | 50 d / 20 ms | 45 d / 97 ms | 50 d / 50 ms | 49 d / 200 ms | 49 d / 244 ms | OK |

### 168-siam100-c3-operator-norm

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c3.machine | 10 | 10 d / 33105 ms | 10 d / 46260 ms | 10 d / 28910 ms | — | — | OK |
| c3.wp50 | 38 | — | — | — | — | — | OK(nobody) |
| c3.machine-n500 | 10 | 15 d / 53 ms | 15 d / 717 ms | 15 d / 49 ms | 15 d / 3080 ms | — | OK |
| c3.wp50-n100 | 48 | 50 d / 234 ms | 50 d / 578 ms | 50 d / 7163 ms | — | — | OK |

### 169-siam100-c4-global-minimum

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c4.machine | 10 | 15 d / 8.3 ms | 15 d / 240 ms | 15 d / 147 ms | — | — | OK |
| c4.wp50 | 48 | 50 d / 10 ms | ✗0 d / 1763 ms | 50 d / 182 ms | — | — | OK |

### 170-siam100-c5-gamma-best-cubic

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c5.machine | 10 | 14 d / 1696 ms | — | 11 d / 19323 ms | — | — | OK |
| c5.wp50 | 48 | 50 d / 2008 ms | — | — | — | — | OK |

### 171-siam100-c6-biased-flea

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c6.machine | 10 | 13 d / 457 ms | ~6 d / 215 ms | 12 d / 844 ms | 12 d / 39625 ms | — | OK |
| c6.wp50 | 48 | — | ~7 d / 2471 ms | — | — | — | FEWER-DIGITS(0<7) |

### 172-siam100-c7-sparse-inverse

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c7.machine | 10 | 15 d / 131 ms | — | 15 d / 2679 ms | — | — | OK |
| c7.wp50 | 48 | 50 d / 679 ms | — | — | — | — | OK |
| c7.machine n=2000 | 10 | 15 d / 57 ms | — | 15 d / 1069 ms | — | 14 d / 160771 ms | OK |
| c7.wp50 n=2000 | 48 | 50 d / 148 ms | — | — | — | — | OK |

### 173-siam100-c8-heat-plate

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c8.machine | 10 | 15 d / 0.27 ms | 15 d / 50 ms | 15 d / 5.8 ms | 15 d / 13 ms | — | OK |
| c8.wp50 | 48 | 50 d / 1.1 ms | 50 d / 3755 ms | 50 d / 10 ms | 50 d / 98 ms | — | OK |
| c8.machine pde | 4 | 4 d / 697 ms | ✗0 d / 760 ms | — | — | — | OK |

### 174-siam100-c9-parametric-maximum

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c9.machine | 10 | 12 d / 510 ms | ~6 d / 98539 ms | ~3 d / 35126 ms | — | — | OK |
| c9.wp50 | 48 | 50 d / 12383 ms | ~30 d / 9224 ms | — | — | — | OK |

### 175-siam100-c10-brownian-ends

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| c10.machine | 10 | 15 d / 0.022 ms | — | 15 d / 0.55 ms | 15 d / 0.85 ms | 15 d / 25 ms | OK |
| c10.wp50 | 48 | 50 d / 0.19 ms | — | 50 d / 1.1 ms | 50 d / 4 ms | 50 d / 25 ms | OK |
| c10.machine pde | 3 | 6 d / 16 ms | — | — | — | — | OK |

### 176-siam100-d1-kempner-42

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d1.machine | 10 | ✗0 d / 799 ms | ✗0 d / 606 ms | ✗0 d / 72 ms | — | — | OK(nobody) |
| d1.wp50 | 48 | — | — | ✗0 d / 1660 ms | — | — | OK(nobody) |

### 177-siam100-d2-digit-length-series

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d2.machine | 10 | ~1 d / 1.6 ms | ~1 d / 27 ms | ~1 d / 59 ms | — | — | OK |
| d2.wp50 | 48 | ~1 d / 2.5 ms | ~1 d / 31 ms | ~1 d / 1614 ms | — | — | OK |

### 178-siam100-d3-thomson-100

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d3.machine | 10 | 14 d / 7512 ms | — | ~9 d / 45152 ms | — | — | OK |
| d3.wp50 | 48 | 50 d / 83087 ms | — | — | — | — | OK |

### 179-siam100-d4-riemann-r-zero

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d4.machine | 10 | ✗0 d / 1.3 ms | — | ✗0 d / 64 ms | — | — | OK(nobody) |
| d4.wp50 | 48 | ✗0 d / 16197 ms | — | ✗0 d / 112 ms | — | — | OK(nobody) |
| d4.wp50-builtin | 48 | 50 d / 14409 ms | — | — | — | — | OK |

### 180-siam100-d5-toeplitz-polynomial

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d5.machine | 10 | ✗0 d / 2099 ms | ✗0 d / 7847 ms | ✗0 d / 866 ms | — | — | OK(nobody) |
| d5.wp50 | 48 | ✗0 d / 347713 ms | ✗0 d / 581693 ms | — | — | — | OK(nobody) |

### 181-siam100-d6-quadruple-exponential

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d6.machine | 10 | 15 d / 0.16 ms | 12 d / 8.5 ms | 15 d / 0.54 ms | 15 d / 11 ms | 15 d / 105 ms | OK |
| d6.wp50 | 48 | 50 d / 75 ms | 50 d / 388 ms | 50 d / 265 ms | ~8 d / 688 ms | 50 d / 152660 ms | OK |

### 182-siam100-d7-region-area

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d7.machine | 10 | 15 d / 6.9 ms | 14 d / 5034 ms | 14 d / 9.7 ms | 15 d / 212 ms | — | OK |
| d7.wp50 | 48 | 50 d / 101 ms | 50 d / 5812 ms | 50 d / 180 ms | ~4 d / 2032 ms | — | OK |

### 183-siam100-d8-sobolev-constant

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d8.machine | 3 | — | — | — | — | — | OK(nobody) |
| d8.wp50 | 3 | — | — | — | — | — | OK(nobody) |
| d8.machine c3R1 | 10 | 15 d / 0.065 ms | 12 d / 3.6 ms | 15 d / 0.14 ms | 15 d / 7 ms | — | OK |
| d8.machine c2R2 | 10 | 15 d / 71 ms | ~9 d / 39 ms | 15 d / 84 ms | 15 d / 1692 ms | — | OK |

### 184-siam100-d9-triangular-walk

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d9.machine | 10 | 14 d / 7.8 ms | 14 d / 26 ms | 14 d / 15 ms | — | 15 d / 50158 ms | FEWER-DIGITS(14<15) |
| d9.wp50 | 48 | 50 d / 562 ms | 50 d / 1831 ms | 50 d / 13765 ms | 50 d / 3509 ms | 50 d / 2103 ms | OK |

### 185-siam100-d10-random-recurrence

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d10.machine | 2 | 4 d / 24 ms | 4 d / 187 ms | 4 d / 242 ms | 4 d / 11807 ms | 3 d / 468 ms | OK |
| d10.wp50 | 1 | 2 d / 24 ms | 2 d / 46 ms | 3 d / 11 ms | 2 d / 1812 ms | 3 d / 35 ms | FEWER-DIGITS(2<3) SLOWER(2.2x) |

### 186-siam100-d11-six-masses

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d11.machine | 10 | 13 d / 55 ms | 13 d / 283 ms | 14 d / 1024 ms | — | — | FEWER-DIGITS(13<14) |
| d11.wp50 | 48 | 50 d / 20686 ms | — | — | — | — | OK |

### 187-siam100-d12-constrained-particle

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d12.machine | 10 | 13 d / 93 ms | 10 d / 276 ms | 12 d / 1899 ms | 11 d / 2409 ms | — | OK |
| d12.wp50 | 48 | 50 d / 1223 ms | — | — | — | — | OK |

### 188-siam100-d13-nonlinear-ode-z0

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d13.machine | 10 | 12 d / 28 ms | 10 d / 227 ms | 11 d / 400 ms | — | — | OK |
| d13.wp50 | 48 | 50 d / 1011 ms | — | 50 d / 178353 ms | — | — | OK |

### 189-siam100-d14-pentagon-poisson

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d14.machine | 10 | ~8 d / 10791 ms | ~3 d / 16308 ms | — | — | — | OK |
| d14.wp50 | 10 | ~8 d / 15011 ms | — | — | — | — | OK(nobody) |

### 190-siam100-d15-blowup-time

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d15.machine | 6 | ~3 d / 3501 ms | ~2 d / 249 ms | ~4 d / 4811 ms | — | — | FEWER-DIGITS(3<4) |
| d15.wp50 | 6 | ~3 d / 3209 ms | ~3 d / 16628 ms | — | — | — | OK |

### 191-siam100-d16-daubechies-heat

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d16.machine | 10 | — | — | — | — | — | OK(nobody) |
| d16.machine-reduced | 10 | ~3 d / 729 ms | — | — | — | — | OK(nobody) |
| d16.wp50 | 48 | — | — | — | — | — | OK(nobody) |

### 192-siam100-d17-lshape-nodal-line

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d17.machine | 10 | ~7 d / 1243 ms | — | — | — | — | OK(nobody) |
| d17.machine-reduced | 10 | ~9 d / 1228 ms | — | — | — | — | OK(nobody) |
| d17.wp50 | 48 | ~7 d / 1242 ms | — | — | — | — | OK(nobody) |

### 193-siam100-d18-koch-eigenvalue

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d18.machine | 5 | ~3 d / 1610 ms | — | — | — | — | OK(nobody) |
| d18.machine level=0 | 10 | ~9 d / 321 ms | — | — | — | — | OK(nobody) |
| d18.wp50 | 5 | ~3 d / 1593 ms | — | — | — | — | OK(nobody) |

### 194-siam100-d19-optimal-design

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d19.machine | 10 | 13 d / 17316 ms | — | ~3 d / 28742 ms | — | — | OK |
| d19.machine-reduced | 10 | 14 d / 166 ms | — | 11 d / 209 ms | 14 d / 5603 ms | — | OK |
| d19.wp50 | 48 | 50 d / 502 ms | — | — | — | — | OK |

### 195-siam100-d20-half-exponential

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| d20.machine | 5 | — | — | — | — | — | OK(nobody) |
| d20.wp50 | 5 | — | — | — | — | — | OK(nobody) |
| d20.machine robbins n=20 | 10 | 14 d / 1.3 ms | — | 15 d / 0.5 ms | 15 d / 19 ms | 15 d / 440 ms | FEWER-DIGITS(14<15) SLOWER(2.6x) |
| d20.wp50 robbins n=60 | 48 | 50 d / 74 ms | — | 48 d / 852 ms | 50 d / 745 ms | 50 d / 10368 ms | OK |

### 196-siam100h-c1-oscillatory-tails

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hc1.1 machine | 10 | 15 d / 0.29 ms | 11 d / 96 ms | ~7 d / 9.5 ms | ~7 d / 78 ms | ~5 d / 183 ms | OK |
| hc1.1 wp50 | 48 | 50 d / 19 ms | ~41 d / 370 ms | ~2 d / 99 ms | — | ~3 d / 2911 ms | OK |
| hc1.2 machine | 10 | 15 d / 0.11 ms | 11 d / 21 ms | ✗0 d / 11 ms | ✗0 d / 103 ms | ✗0 d / 159 ms | OK |
| hc1.2 wp50 | 48 | 50 d / 19 ms | ~43 d / 289 ms | ✗0 d / 56 ms | — | ✗0 d / 3085 ms | OK |
| hc1.3 machine | 10 | 15 d / 0.13 ms | ~9 d / 18 ms | ✗0 d / 12 ms | ✗0 d / 104 ms | ✗0 d / 161 ms | OK |
| hc1.3 wp50 | 48 | 50 d / 19 ms | 49 d / 714 ms | ✗0 d / 52 ms | — | ✗0 d / 3312 ms | OK |
| hc1.4 machine | 10 | 15 d / 0.14 ms | ~7 d / 42 ms | ~2 d / 12 ms | ~2 d / 126 ms | ✗0 d / 330 ms | OK |
| hc1.4 wp50 | 48 | 50 d / 21 ms | ~7 d / 267 ms | ~1 d / 55 ms | — | ✗0 d / 3576 ms | OK |
| hc1.5 machine | 10 | 12 d / 0.34 ms | ~7 d / 100 ms | ✗0 d / 14 ms | ✗0 d / 171 ms | ✗0 d / 450 ms | OK |
| hc1.5 wp50 | 48 | 50 d / 159 ms | ~9 d / 295 ms | ✗0 d / 72 ms | — | ✗0 d / 6520 ms | OK |
| hc1.6 machine | 10 | 12 d / 0.21 ms | ~8 d / 23 ms | ~3 d / 13 ms | ~3 d / 155 ms | ~4 d / 310 ms | OK |
| hc1.6 wp50 | 48 | 50 d / 146 ms | ~13 d / 281 ms | ~2 d / 62 ms | — | ~4 d / 6369 ms | OK |
| hc1.7 machine | 10 | 12 d / 3.1 ms | ✗0 d / 119 ms | — | — | ~3 d / 178 ms | OK |
| hc1.7 wp50 | 48 | 50 d / 376 ms | ~1 d / 468 ms | — | — | ~1 d / 6406 ms | OK |
| hc1.8 machine | 10 | 12 d / 0.38 ms | 12 d / 61 ms | ~5 d / 13 ms | ~5 d / 124 ms | ~4 d / 126 ms | OK |
| hc1.8 wp50 | 48 | 50 d / 292 ms | ~44 d / 884 ms | ~2 d / 42 ms | — | ~3 d / 3503 ms | OK |
| hc1.9 machine | 10 | 12 d / 0.2 ms | 11 d / 27 ms | ✗0 d / 10 ms | ✗0 d / 143 ms | ✗0 d / 110 ms | OK |
| hc1.9 wp50 | 48 | 50 d / 173 ms | ~40 d / 275 ms | ✗0 d / 48 ms | — | ✗0 d / 3378 ms | OK |
| hc1.10 machine | 10 | 12 d / 0.29 ms | ~3 d / 76 ms | ✗0 d / 12 ms | ✗0 d / 195 ms | ✗0 d / 464 ms | OK |
| hc1.10 wp50 | 48 | 50 d / 157 ms | ~3 d / 299 ms | ✗0 d / 61 ms | — | ✗0 d / 7098 ms | OK |

### 197-siam100h-c2-mirror-billiards

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hc2.1 machine | 10 | 15 d / 0.86 ms | 15 d / 78 ms | 15 d / 9.1 ms | 15 d / 1.4 ms | 15 d / 119 ms | OK |
| hc2.1 wp50 | 48 | 50 d / 1.3 ms | 50 d / 7.8 ms | 50 d / 2.8 ms | 50 d / 8.5 ms | 50 d / 106 ms | OK |
| hc2.2 machine | 10 | 12 d / 1.6 ms | 12 d / 10 ms | 12 d / 1.3 ms | 13 d / 3.1 ms | 12 d / 10 ms | FEWER-DIGITS(12<13) |
| hc2.2 wp50 | 48 | 50 d / 4.6 ms | 50 d / 21 ms | 50 d / 13 ms | 50 d / 34 ms | 50 d / 39 ms | OK |
| hc2.3 machine | 10 | 13 d / 1.3 ms | 13 d / 7.1 ms | 14 d / 0.72 ms | 14 d / 2.6 ms | 13 d / 7.8 ms | FEWER-DIGITS(13<14) SLOWER(1.8x) |
| hc2.3 wp50 | 48 | 50 d / 3.7 ms | 50 d / 16 ms | 50 d / 7.1 ms | 50 d / 30 ms | 50 d / 25 ms | OK |
| hc2.4 machine | 10 | 12 d / 2.8 ms | 12 d / 17 ms | 11 d / 1.9 ms | 11 d / 5.5 ms | 11 d / 18 ms | OK |
| hc2.4 wp50 | 48 | 50 d / 7.3 ms | 50 d / 34 ms | 50 d / 17 ms | 50 d / 81 ms | 50 d / 65 ms | OK |
| hc2.5 machine | 10 | 15 d / 0.83 ms | 14 d / 4.8 ms | 14 d / 0.49 ms | 15 d / 1.5 ms | 15 d / 4.9 ms | SLOWER(1.7x) |
| hc2.5 wp50 | 48 | 50 d / 2.1 ms | 50 d / 10 ms | 50 d / 5.1 ms | 50 d / 13 ms | 50 d / 20 ms | OK |
| hc2.6 machine | 10 | 12 d / 2.5 ms | 12 d / 13 ms | 12 d / 4 ms | 12 d / 5 ms | 12 d / 14 ms | OK |
| hc2.6 wp50 | 48 | 50 d / 5.6 ms | 50 d / 24 ms | 50 d / 14 ms | 50 d / 54 ms | 50 d / 50 ms | OK |
| hc2.7 machine | 6 | 7 d / 3.2 ms | 7 d / 18 ms | 7 d / 1.6 ms | 7 d / 7.2 ms | 7 d / 16 ms | SLOWER(2.0x) |
| hc2.7 wp50 | 48 | 50 d / 11 ms | 50 d / 45 ms | 50 d / 18 ms | 50 d / 62 ms | 50 d / 56 ms | OK |
| hc2.8 machine | 10 | 13 d / 1.3 ms | 13 d / 9.2 ms | 13 d / 0.8 ms | 13 d / 2.4 ms | 13 d / 8.9 ms | SLOWER(1.7x) |
| hc2.8 wp50 | 48 | 50 d / 3.5 ms | 50 d / 15 ms | 50 d / 12 ms | 50 d / 27 ms | 50 d / 29 ms | OK |
| hc2.9 machine | 10 | 13 d / 1.1 ms | 12 d / 4.2 ms | 13 d / 0.59 ms | 12 d / 2.8 ms | 13 d / 4 ms | SLOWER(1.9x) |
| hc2.9 wp50 | 48 | 50 d / 2.6 ms | 50 d / 10 ms | 50 d / 4.7 ms | 50 d / 15 ms | 50 d / 14 ms | OK |
| hc2.10 machine | 10 | 12 d / 2.2 ms | 11 d / 14 ms | 11 d / 1.3 ms | 12 d / 4.4 ms | 12 d / 15 ms | SLOWER(1.7x) |
| hc2.10 wp50 | 48 | 50 d / 6.3 ms | 50 d / 28 ms | 50 d / 15 ms | 50 d / 49 ms | 50 d / 47 ms | OK |

### 198-siam100h-c3-matrix-norms

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hc3.1 machine | 10 | 15 d / 11 ms | 15 d / 135 ms | 15 d / 5.4 ms | 15 d / 54 ms | — | SLOWER(2.0x) |
| hc3.1 wp50 | 48 | 50 d / 103 ms | 50 d / 544 ms | 50 d / 4133 ms | — | — | OK |
| hc3.2 machine | 10 | 15 d / 0.52 ms | 15 d / 2.8 ms | 15 d / 0.48 ms | 15 d / 5.3 ms | — | OK |
| hc3.2 wp50 | 48 | 50 d / 4.9 ms | 50 d / 16 ms | 50 d / 117 ms | — | — | OK |
| hc3.3 machine | 10 | 15 d / 1.2 ms | 15 d / 3.8 ms | 15 d / 0.92 ms | 15 d / 9.6 ms | — | OK |
| hc3.3 wp50 | 48 | 50 d / 20 ms | 50 d / 33 ms | 50 d / 277 ms | — | — | OK |
| hc3.4 machine | 10 | 15 d / 9.1 ms | 15 d / 119 ms | 15 d / 2.8 ms | 15 d / 54 ms | — | SLOWER(3.3x) |
| hc3.4 wp50 | 48 | 50 d / 222 ms | 50 d / 544 ms | 50 d / 4182 ms | — | — | OK |
| hc3.5 machine | 10 | 15 d / 0.71 ms | 15 d / 23 ms | 15 d / 2.2 ms | 15 d / 44 ms | — | OK |
| hc3.5 wp50 | 48 | 50 d / 221 ms | 50 d / 466 ms | 50 d / 4368 ms | — | — | OK |
| hc3.6 machine | 10 | 15 d / 0.94 ms | 15 d / 36 ms | 15 d / 4.6 ms | 15 d / 93 ms | — | OK |
| hc3.6 wp50 | 48 | 50 d / 119 ms | 50 d / 540 ms | 50 d / 4736 ms | — | — | OK |
| hc3.7 machine | 10 | 15 d / 1.2 ms | 15 d / 48 ms | 15 d / 3.8 ms | 15 d / 74 ms | — | OK |
| hc3.7 wp50 | 48 | 50 d / 231 ms | 50 d / 499 ms | 50 d / 4440 ms | — | — | OK |
| hc3.8 machine | 10 | 15 d / 0.53 ms | 15 d / 21 ms | 15 d / 1.8 ms | 15 d / 31 ms | — | OK |
| hc3.8 wp50 | 48 | 50 d / 117 ms | 50 d / 244 ms | 50 d / 2401 ms | — | — | OK |
| hc3.9 machine | 10 | 15 d / 2 ms | 15 d / 83 ms | 15 d / 7.4 ms | 15 d / 165 ms | — | OK |
| hc3.9 wp50 | 48 | 50 d / 341 ms | 50 d / 1497 ms | 50 d / 15795 ms | — | — | OK |
| hc3.10 machine | 10 | 15 d / 1.6 ms | 15 d / 40 ms | 15 d / 3.2 ms | 15 d / 84 ms | — | OK |
| hc3.10 wp50 | 48 | 50 d / 167 ms | 50 d / 562 ms | 50 d / 5160 ms | — | — | OK |

### 199-siam100h-c4-global-minima

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hc4.1 machine | 10 | 15 d / 1.3 ms | 15 d / 436 ms | 13 d / 66 ms | — | — | OK |
| hc4.1 wp50 | 48 | 50 d / 2 ms | ✗0 d / 174 ms | 50 d / 58 ms | — | — | OK |
| hc4.2 machine | 10 | 15 d / 3.4 ms | 15 d / 19 ms | 13 d / 59 ms | — | — | OK |
| hc4.2 wp50 | 48 | 50 d / 5 ms | — | 50 d / 62 ms | — | — | OK |
| hc4.3 machine | 10 | 15 d / 1.5 ms | 15 d / 2258 ms | 15 d / 112 ms | — | — | OK |
| hc4.3 wp50 | 48 | 50 d / 1.8 ms | — | 50 d / 122 ms | — | — | OK |
| hc4.4 machine | 10 | 15 d / 12 ms | 15 d / 3728 ms | 15 d / 111 ms | — | — | OK |
| hc4.4 wp50 | 48 | 50 d / 15 ms | — | 50 d / 134 ms | — | — | OK |
| hc4.5 machine | 10 | 15 d / 1.7 ms | 15 d / 2272 ms | 13 d / 63 ms | — | — | OK |
| hc4.5 wp50 | 48 | 50 d / 4 ms | — | 50 d / 70 ms | — | — | OK |
| hc4.6 machine | 10 | 15 d / 4.2 ms | 15 d / 18 ms | 14 d / 60 ms | — | — | OK |
| hc4.6 wp50 | 48 | 50 d / 5.3 ms | ✗0 d / 589 ms | 50 d / 71 ms | — | — | OK |
| hc4.7 machine | 10 | 15 d / 7.9 ms | 15 d / 2367 ms | ~1 d / 84 ms | — | — | OK |
| hc4.7 wp50 | 48 | 50 d / 9.4 ms | — | ~1 d / 97 ms | — | — | OK |
| hc4.8 machine | 10 | 15 d / 2.7 ms | ~7 d / 127 ms | 13 d / 131 ms | — | — | OK |
| hc4.8 wp50 | 48 | ~20 d / 4.4 ms | ✗0 d / 66 ms | ~1 d / 142 ms | — | — | OK |
| hc4.9 machine | 10 | 15 d / 0.94 ms | ✗0 d / 123 ms | ✗0 d / 71 ms | — | — | OK |
| hc4.9 wp50 | 48 | 50 d / 1.4 ms | ✗0 d / 48 ms | ✗0 d / 74 ms | — | — | OK |
| hc4.10 machine | 10 | 15 d / 10 ms | ✗0 d / 1488 ms | ✗0 d / 97 ms | — | — | OK |
| hc4.10 wp50 | 48 | 50 d / 12 ms | ✗0 d / 842 ms | ✗0 d / 118 ms | — | — | OK |

### 200-siam100h-c5-discrete-minimax

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hc5.1 machine | 10 | 15 d / 23 ms | — | 12 d / 215 ms | — | — | OK |
| hc5.1 wp50 | 48 | 50 d / 27 ms | — | — | — | — | OK |
| hc5.2 machine | 10 | 13 d / 52 ms | — | 12 d / 1127 ms | — | — | OK |
| hc5.2 wp50 | 48 | 50 d / 58 ms | — | — | — | — | OK |
| hc5.3 machine | 10 | 15 d / 197 ms | — | 12 d / 2168 ms | — | — | OK |
| hc5.3 wp50 | 48 | 50 d / 211 ms | ✗0 d / 13297 ms | — | — | — | OK |
| hc5.4 machine | 10 | 13 d / 86 ms | — | 12 d / 1517 ms | — | — | OK |
| hc5.4 wp50 | 48 | 50 d / 94 ms | ✗0 d / 10921 ms | — | — | — | OK |
| hc5.5 machine | 10 | 14 d / 285 ms | — | 12 d / 4554 ms | — | — | OK |
| hc5.5 wp50 | 48 | 50 d / 301 ms | — | — | — | — | OK |
| hc5.6 machine | 10 | 14 d / 26 ms | 14 d / 12 ms | 12 d / 564 ms | — | — | SLOWER(2.2x) |
| hc5.6 wp50 | 48 | 50 d / 30 ms | ✗0 d / 1486 ms | — | — | — | OK |
| hc5.7 machine | 10 | 15 d / 389 ms | ~6 d / 6710 ms | 12 d / 2941 ms | — | — | OK |
| hc5.7 wp50 | 48 | 50 d / 415 ms | ~2 d / 47061 ms | — | — | — | OK |
| hc5.8 machine | 10 | 14 d / 652 ms | ~4 d / 7061 ms | 11 d / 3929 ms | — | — | OK |
| hc5.8 wp50 | 48 | 50 d / 626 ms | ~3 d / 45324 ms | — | — | — | OK |
| hc5.9 machine | 10 | 14 d / 194 ms | ~5 d / 4207 ms | 12 d / 3796 ms | — | — | OK |
| hc5.9 wp50 | 48 | 50 d / 220 ms | ~3 d / 36235 ms | — | — | — | OK |
| hc5.10 machine | 10 | 15 d / 2292 ms | ~5 d / 7431 ms | 12 d / 26461 ms | — | — | OK |
| hc5.10 wp50 | 48 | 50 d / 2401 ms | — | — | — | — | OK |

### 201-siam100h-c6-lattice-return

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hc6.1 machine | 10 | 13 d / 32 ms | ~8 d / 127 ms | 13 d / 248 ms | 13 d / 5293 ms | — | OK |
| hc6.1 wp50 | 48 | — | ~8 d / 341 ms | — | — | — | FEWER-DIGITS(0<8) |
| hc6.2 machine | 10 | 13 d / 240 ms | ~8 d / 199 ms | 13 d / 1846 ms | — | — | OK |
| hc6.2 wp50 | 48 | — | ~9 d / 1991 ms | — | — | — | FEWER-DIGITS(0<9) |
| hc6.3 machine | 10 | 14 d / 397 ms | ~7 d / 363 ms | 12 d / 5104 ms | — | — | OK |
| hc6.3 wp50 | 48 | — | ~7 d / 4446 ms | — | — | — | FEWER-DIGITS(0<7) |
| hc6.4 machine | 10 | 11 d / 1858 ms | ~2 d / 465 ms | — | — | — | OK |
| hc6.4 wp50 | 48 | — | ~2 d / 45832 ms | — | — | — | FEWER-DIGITS(0<2) |
| hc6.5 machine | 10 | 13 d / 456 ms | 10 d / 13059 ms | 11 d / 3426 ms | — | — | OK |
| hc6.5 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hc6.6 machine | 10 | 13 d / 228 ms | ~7 d / 6180 ms | 13 d / 1062 ms | 13 d / 41369 ms | — | OK |
| hc6.6 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hc6.7 machine | 10 | 13 d / 1703 ms | ~6 d / 4961 ms | 12 d / 21031 ms | — | — | OK |
| hc6.7 wp50 | 48 | — | ~6 d / 38197 ms | — | — | — | FEWER-DIGITS(0<6) |
| hc6.8 machine | 10 | 13 d / 4604 ms | ~8 d / 22208 ms | — | — | — | OK |
| hc6.8 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hc6.9 machine | 10 | 14 d / 184 ms | ~2 d / 40 ms | ~3 d / 6224 ms | — | — | OK |
| hc6.9 wp50 | 48 | — | ~2 d / 290 ms | — | — | — | FEWER-DIGITS(0<2) |
| hc6.10 machine | 10 | 13 d / 1263 ms | ~5 d / 179 ms | 12 d / 14208 ms | — | — | OK |
| hc6.10 wp50 | 48 | — | ~5 d / 2327 ms | — | — | — | FEWER-DIGITS(0<5) |

### 202-siam100h-c7-sparse-solve

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hc7.1 machine | 10 | 15 d / 1.1 ms | 15 d / 45 ms | 15 d / 37 ms | — | — | OK |
| hc7.1 wp50 | 48 | 50 d / 28 ms | 50 d / 1147 ms | — | — | — | OK |
| hc7.2 machine | 10 | 15 d / 0.76 ms | 15 d / 10 ms | 15 d / 20 ms | — | — | OK |
| hc7.2 wp50 | 48 | 50 d / 14 ms | 50 d / 565 ms | — | — | — | OK |
| hc7.3 machine | 10 | 15 d / 138 ms | 15 d / 6672 ms | — | — | — | OK |
| hc7.3 wp50 | 48 | 50 d / 217 ms | 50 d / 5545 ms | — | — | — | OK |
| hc7.4 machine | 10 | 15 d / 242 ms | 15 d / 5589 ms | 15 d / 10045 ms | — | — | OK |
| hc7.4 wp50 | 48 | 50 d / 478 ms | 50 d / 3763 ms | — | — | — | OK |
| hc7.5 machine | 10 | 15 d / 9 ms | 14 d / 247 ms | 14 d / 253 ms | 15 d / 46076 ms | 14 d / 17137 ms | OK |
| hc7.5 wp50 | 48 | 50 d / 32 ms | 50 d / 338 ms | — | — | — | OK |
| hc7.6 machine | 10 | 15 d / 228 ms | 15 d / 2537 ms | 15 d / 7346 ms | — | — | OK |
| hc7.6 wp50 | 48 | 50 d / 457 ms | 50 d / 4193 ms | — | — | — | OK |
| hc7.7 machine | 10 | 15 d / 9 ms | 15 d / 111 ms | 15 d / 191 ms | 15 d / 39853 ms | 15 d / 20917 ms | OK |
| hc7.7 wp50 | 48 | 50 d / 31 ms | 50 d / 1279 ms | — | — | — | OK |
| hc7.8 machine | 10 | 15 d / 1.7 ms | 15 d / 46 ms | 15 d / 12 ms | — | — | OK |
| hc7.8 wp50 | 48 | 50 d / 16 ms | 50 d / 791 ms | — | — | — | OK |
| hc7.9 machine | 10 | 15 d / 10 ms | 15 d / 111 ms | 15 d / 158 ms | 15 d / 37141 ms | 15 d / 19242 ms | OK |
| hc7.9 wp50 | 48 | 50 d / 26 ms | 50 d / 780 ms | — | — | — | OK |
| hc7.10 machine | 10 | 15 d / 57 ms | 15 d / 619 ms | 15 d / 919 ms | — | — | OK |
| hc7.10 wp50 | 48 | 50 d / 109 ms | 50 d / 3049 ms | — | — | — | OK |

### 203-siam100h-c8-heat-crossing

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hc8.1 machine | 10 | 15 d / 0.27 ms | 15 d / 57 ms | 15 d / 7.2 ms | 15 d / 10 ms | — | OK |
| hc8.1 wp50 | 48 | 50 d / 1.2 ms | 50 d / 4605 ms | 50 d / 14 ms | 50 d / 231 ms | — | OK |
| hc8.2 machine | 10 | 15 d / 0.37 ms | 15 d / 67 ms | 15 d / 7.1 ms | 15 d / 7.8 ms | — | OK |
| hc8.2 wp50 | 48 | 50 d / 1.2 ms | 50 d / 5285 ms | 50 d / 10 ms | 50 d / 283 ms | — | OK |
| hc8.3 machine | 10 | 15 d / 0.25 ms | 15 d / 70 ms | 15 d / 6.8 ms | 15 d / 9.4 ms | — | OK |
| hc8.3 wp50 | 48 | 50 d / 0.95 ms | 50 d / 4388 ms | 50 d / 14 ms | 50 d / 299 ms | — | OK |
| hc8.4 machine | 10 | 15 d / 0.23 ms | 15 d / 32 ms | 15 d / 5.2 ms | 15 d / 9.6 ms | — | OK |
| hc8.4 wp50 | 48 | 50 d / 0.88 ms | 50 d / 393 ms | 50 d / 9.4 ms | 50 d / 245 ms | — | OK |
| hc8.5 machine | 10 | 15 d / 0.31 ms | 15 d / 70 ms | 14 d / 6.9 ms | 15 d / 9.8 ms | — | OK |
| hc8.5 wp50 | 48 | 50 d / 1.2 ms | 50 d / 5126 ms | 50 d / 14 ms | 50 d / 309 ms | — | OK |
| hc8.6 machine | 10 | 15 d / 0.58 ms | 15 d / 162 ms | 15 d / 17 ms | 15 d / 14 ms | — | OK |
| hc8.6 wp50 | 48 | 50 d / 2.3 ms | 50 d / 10924 ms | 50 d / 50 ms | 50 d / 421 ms | — | OK |
| hc8.7 machine | 10 | 15 d / 0.91 ms | — | 15 d / 147 ms | — | — | OK |
| hc8.7 wp50 | 48 | 50 d / 128 ms | — | 50 d / 642 ms | — | — | OK |
| hc8.8 machine | 10 | 15 d / 0.18 ms | — | 15 d / 6.5 ms | 15 d / 7.5 ms | — | OK |
| hc8.8 wp50 | 48 | 50 d / 1 ms | — | 50 d / 39 ms | 50 d / 269 ms | — | OK |
| hc8.9 machine | 10 | 15 d / 0.62 ms | — | 15 d / 12 ms | 15 d / 36 ms | — | OK |
| hc8.9 wp50 | 48 | 50 d / 2.9 ms | — | 50 d / 84 ms | 50 d / 730 ms | — | OK |
| hc8.10 machine | 10 | 14 d / 0.34 ms | — | 14 d / 14 ms | 14 d / 9.4 ms | — | OK |
| hc8.10 wp50 | 48 | 50 d / 2 ms | — | 50 d / 36 ms | 50 d / 308 ms | — | OK |

### 204-siam100h-c9-parametric-extremum

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hc9.1 machine | 10 | 12 d / 477 ms | — | ~2 d / 12080 ms | — | — | OK |
| hc9.1 wp50 | 48 | 50 d / 11856 ms | — | — | — | — | OK |
| hc9.2 machine | 10 | 10 d / 563 ms | — | ~2 d / 13403 ms | — | — | OK |
| hc9.2 wp50 | 48 | ~33 d / 8854 ms | — | — | — | — | OK(nobody) |
| hc9.3 machine | 10 | 13 d / 46 ms | — | ~6 d / 38 ms | — | — | OK |
| hc9.3 wp50 | 48 | 50 d / 848 ms | — | 50 d / 2511 ms | — | — | OK |
| hc9.4 machine | 10 | 13 d / 161 ms | — | ~8 d / 24 ms | — | — | OK |
| hc9.4 wp50 | 48 | 50 d / 571 ms | — | 50 d / 943 ms | — | — | OK |
| hc9.5 machine | 10 | 10 d / 131 ms | — | ~6 d / 61 ms | — | — | OK |
| hc9.5 wp50 | 48 | 50 d / 338 ms | — | 50 d / 10030 ms | — | — | OK |
| hc9.6 machine | 10 | 11 d / 43 ms | — | ~7 d / 20 ms | — | — | OK |
| hc9.6 wp50 | 48 | 50 d / 146 ms | — | 50 d / 303 ms | — | — | OK |
| hc9.7 machine | 10 | 12 d / 154 ms | — | ~7 d / 52 ms | — | — | OK |
| hc9.7 wp50 | 48 | 50 d / 20544 ms | — | 50 d / 8070 ms | — | — | SLOWER(2.5x) |
| hc9.8 machine | 10 | 12 d / 166 ms | ~8 d / 22680 ms | ~7 d / 34 ms | — | — | OK |
| hc9.8 wp50 | 48 | 50 d / 1823 ms | — | 50 d / 2170 ms | — | — | OK |
| hc9.9 machine | 10 | 10 d / 1350 ms | — | ~3 d / 12982 ms | — | — | OK |
| hc9.9 wp50 | 48 | 48 d / 10956 ms | ~23 d / 6476 ms | — | — | — | OK |
| hc9.10 machine | 10 | 12 d / 507 ms | ~9 d / 59947 ms | ~3 d / 8115 ms | — | — | OK |
| hc9.10 wp50 | 48 | 50 d / 11702 ms | ~30 d / 6279 ms | — | — | — | OK |

### 205-siam100h-c10-brownian-hitting

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hc10.1 machine | 10 | 15 d / 0.025 ms | 15 d / 58 ms | 15 d / 0.39 ms | 15 d / 0.17 ms | 15 d / 31 ms | OK |
| hc10.1 wp50 | 48 | 50 d / 0.23 ms | 50 d / 61 ms | 50 d / 2.1 ms | 50 d / 2 ms | 50 d / 18 ms | OK |
| hc10.2 machine | 10 | 15 d / 0.098 ms | 15 d / 179 ms | 15 d / 0.52 ms | 15 d / 0.21 ms | 15 d / 8.1 ms | OK |
| hc10.2 wp50 | 48 | 50 d / 0.4 ms | 50 d / 183 ms | 50 d / 2.6 ms | 50 d / 2.8 ms | 50 d / 17 ms | OK |
| hc10.3 machine | 10 | 15 d / 0.17 ms | — | 15 d / 0.72 ms | 15 d / 0.4 ms | 15 d / 4.8 ms | OK |
| hc10.3 wp50 | 48 | 50 d / 0.81 ms | ~42 d / 1931 ms | 50 d / 3.4 ms | 50 d / 5.6 ms | 50 d / 10 ms | OK |
| hc10.4 machine | 10 | 15 d / 0.027 ms | 15 d / 72 ms | 15 d / 1.3 ms | 15 d / 0.66 ms | 15 d / 4.2 ms | OK |
| hc10.4 wp50 | 48 | 50 d / 0.47 ms | ~30 d / 73 ms | 50 d / 11 ms | 50 d / 19 ms | 50 d / 8.9 ms | OK |
| hc10.5 machine | 10 | 15 d / 0.037 ms | 15 d / 711 ms | 15 d / 2.8 ms | 15 d / 1.4 ms | 15 d / 8.9 ms | OK |
| hc10.5 wp50 | 48 | 50 d / 3.7 ms | ~31 d / 717 ms | 50 d / 18 ms | 50 d / 24 ms | 50 d / 28 ms | OK |
| hc10.6 machine | 10 | 15 d / 0.021 ms | 15 d / 70 ms | 15 d / 0.26 ms | 15 d / 0.14 ms | 15 d / 4.1 ms | OK |
| hc10.6 wp50 | 48 | 50 d / 0.12 ms | 50 d / 70 ms | 50 d / 1.3 ms | 50 d / 1.6 ms | 50 d / 7.1 ms | OK |
| hc10.7 machine | 10 | 15 d / 0.024 ms | 15 d / 74 ms | 15 d / 0.39 ms | 15 d / 0.18 ms | 15 d / 7.1 ms | OK |
| hc10.7 wp50 | 48 | 50 d / 0.14 ms | 50 d / 74 ms | 50 d / 1.7 ms | 50 d / 3 ms | 50 d / 8.8 ms | OK |
| hc10.8 machine | 10 | 15 d / 0.16 ms | — | 15 d / 4.9 ms | 15 d / 33 ms | 15 d / 13 ms | OK |
| hc10.8 wp50 | 48 | 50 d / 2.1 ms | — | 50 d / 54 ms | 50 d / 86 ms | 50 d / 27 ms | OK |
| hc10.9 machine | 10 | 15 d / 0.83 ms | — | 15 d / 55 ms | 15 d / 184 ms | 15 d / 20 ms | OK |
| hc10.9 wp50 | 48 | 50 d / 129 ms | — | 50 d / 326 ms | 50 d / 2072 ms | 50 d / 664 ms | OK |
| hc10.10 machine | 10 | 15 d / 1.6 ms | 15 d / 153 ms | 15 d / 3.4 ms | 15 d / 3.2 ms | 15 d / 6.2 ms | OK |
| hc10.10 wp50 | 48 | 50 d / 6.2 ms | ~21 d / 153 ms | 50 d / 31 ms | 50 d / 50 ms | 50 d / 21 ms | OK |

### 206-siam100h-d1-kempner-family

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd1.1 machine | 10 | ✗0 d / 884 ms | ✗0 d / 379 ms | ✗0 d / 47 ms | — | — | OK(nobody) |
| hd1.1 wp50 | 48 | — | ✗0 d / 57 ms | ✗0 d / 1076 ms | — | — | OK(nobody) |
| hd1.2 machine | 10 | ✗0 d / 829 ms | ✗0 d / 265 ms | ✗0 d / 0.15 ms | — | — | OK(nobody) |
| hd1.2 wp50 | 48 | — | — | ✗0 d / 0.12 ms | — | — | OK(nobody) |
| hd1.3 machine | 10 | ✗0 d / 948 ms | ✗0 d / 359 ms | ✗0 d / 8.1 ms | — | — | OK(nobody) |
| hd1.3 wp50 | 48 | — | — | ✗0 d / 8.9 ms | — | — | OK(nobody) |
| hd1.4 machine | 10 | ✗0 d / 981 ms | ✗0 d / 197 ms | ✗0 d / 47 ms | — | — | OK(nobody) |
| hd1.4 wp50 | 48 | — | — | ✗0 d / 994 ms | — | — | OK(nobody) |
| hd1.5 machine | 10 | ✗0 d / 699 ms | ✗0 d / 251 ms | ✗0 d / 0.44 ms | — | — | OK(nobody) |
| hd1.5 wp50 | 48 | — | — | ✗0 d / 0.4 ms | — | — | OK(nobody) |
| hd1.6 machine | 10 | ✗0 d / 0.25 ms | ✗0 d / 224 ms | ✗0 d / 0.38 ms | — | — | OK(nobody) |
| hd1.6 wp50 | 48 | ✗0 d / 0.56 ms | ✗0 d / 53 ms | ✗0 d / 0.38 ms | — | — | OK(nobody) |
| hd1.7 machine | 10 | ~1 d / 32 ms | ~3 d / 182 ms | ~2 d / 48 ms | — | — | FEWER-DIGITS(1<3) |
| hd1.7 wp50 | 48 | — | ~2 d / 53 ms | ~3 d / 1007 ms | — | — | FEWER-DIGITS(0<3) |
| hd1.8 machine | 10 | ✗0 d / 3.1 ms | ✗0 d / 183 ms | ✗0 d / 0.14 ms | — | — | OK(nobody) |
| hd1.8 wp50 | 48 | ✗0 d / 37975 ms | ✗0 d / 51 ms | ✗0 d / 0.11 ms | — | — | OK(nobody) |
| hd1.9 machine | 10 | ✗0 d / 0.27 ms | ~1 d / 203 ms | ✗0 d / 0.12 ms | — | — | FEWER-DIGITS(0<1) |
| hd1.9 wp50 | 48 | ~2 d / 0.51 ms | ~1 d / 53 ms | ✗0 d / 0.12 ms | — | — | OK |
| hd1.10 machine | 10 | ✗0 d / 1158 ms | ✗0 d / 192 ms | ✗0 d / 44 ms | — | — | OK(nobody) |
| hd1.10 wp50 | 48 | — | — | ✗0 d / 997 ms | — | — | OK(nobody) |

### 207-siam100h-d2-digit-length-family

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd2.1 machine | 10 | ~7 d / 2.4 ms | ~4 d / 293 ms | ~6 d / 50 ms | — | — | OK |
| hd2.1 wp50 | 48 | ~7 d / 3.1 ms | ~4 d / 22 ms | ~7 d / 1047 ms | — | — | OK |
| hd2.2 machine | 10 | ~1 d / 2.1 ms | ~1 d / 201 ms | ~1 d / 44 ms | — | — | OK |
| hd2.2 wp50 | 48 | ~1 d / 3.2 ms | ~1 d / 21 ms | ~1 d / 1002 ms | — | — | OK |
| hd2.3 machine | 10 | ~3 d / 0.45 ms | ~3 d / 404 ms | ~3 d / 47 ms | — | — | OK |
| hd2.3 wp50 | 48 | ~3 d / 0.74 ms | ~3 d / 34 ms | ~5 d / 728 ms | — | — | FEWER-DIGITS(3<5) |
| hd2.4 machine | 10 | ✗0 d / 1.1 ms | ~1 d / 203 ms | ✗0 d / 43 ms | — | — | FEWER-DIGITS(0<1) |
| hd2.4 wp50 | 48 | ✗0 d / 1.6 ms | ~1 d / 22 ms | ~1 d / 1011 ms | — | — | FEWER-DIGITS(0<1) |
| hd2.5 machine | 10 | ~2 d / 1.4 ms | ~1 d / 202 ms | ~2 d / 49 ms | — | — | OK |
| hd2.5 wp50 | 48 | ~2 d / 2.1 ms | ~1 d / 22 ms | ~2 d / 1028 ms | — | — | OK |
| hd2.6 machine | 10 | ~7 d / 1.7 ms | ~4 d / 154 ms | ~6 d / 46 ms | — | — | OK |
| hd2.6 wp50 | 48 | ~7 d / 2.4 ms | ~4 d / 34 ms | ~8 d / 1069 ms | — | — | FEWER-DIGITS(7<8) |
| hd2.7 machine | 10 | ✗0 d / 1.6 ms | ✗0 d / 204 ms | ✗0 d / 46 ms | — | — | OK(nobody) |
| hd2.7 wp50 | 48 | ✗0 d / 2.4 ms | ✗0 d / 21 ms | ✗0 d / 1000 ms | — | — | OK(nobody) |
| hd2.8 machine | 10 | ~3 d / 0.25 ms | ~3 d / 172 ms | ~3 d / 7.7 ms | — | — | OK |
| hd2.8 wp50 | 48 | ~3 d / 0.52 ms | ~3 d / 32 ms | ~3 d / 26 ms | — | — | OK |
| hd2.9 machine | 10 | ~2 d / 6.6 ms | ~1 d / 210 ms | ~2 d / 47 ms | — | — | OK |
| hd2.9 wp50 | 48 | ~2 d / 11 ms | ~1 d / 23 ms | ~2 d / 1003 ms | — | — | OK |
| hd2.10 machine | 10 | ~1 d / 1.8 ms | ~1 d / 220 ms | ~1 d / 48 ms | — | — | OK |
| hd2.10 wp50 | 48 | ~1 d / 2.6 ms | ~1 d / 28 ms | ~1 d / 1020 ms | — | — | OK |

### 208-siam100h-d3-point-energy-family

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd3.1 machine | 10 | 15 d / 4.4 ms | 15 d / 15 ms | 11 d / 20 ms | 15 d / 147 ms | — | OK |
| hd3.1 wp50 | 48 | 50 d / 30 ms | 50 d / 44 ms | — | — | — | OK |
| hd3.2 machine | 10 | 15 d / 24 ms | 15 d / 72 ms | 12 d / 57 ms | 15 d / 1205 ms | — | OK |
| hd3.2 wp50 | 48 | 50 d / 204 ms | 50 d / 242 ms | — | — | — | OK |
| hd3.3 machine | 10 | 15 d / 8 ms | 15 d / 11 ms | 11 d / 22 ms | 15 d / 479 ms | — | OK |
| hd3.3 wp50 | 48 | 50 d / 67 ms | 50 d / 104 ms | — | — | — | OK |
| hd3.4 machine | 10 | 14 d / 355 ms | 10 d / 265 ms | 12 d / 452 ms | 13 d / 35164 ms | — | OK |
| hd3.4 wp50 | 48 | 50 d / 3407 ms | ~10 d / 2317 ms | — | — | — | OK |
| hd3.5 machine | 10 | 15 d / 94 ms | 15 d / 82 ms | 12 d / 225 ms | 15 d / 3896 ms | — | OK |
| hd3.5 wp50 | 48 | 50 d / 1008 ms | ~16 d / 756 ms | — | — | — | OK |
| hd3.6 machine | 10 | 15 d / 328 ms | ~5 d / 218 ms | 11 d / 723 ms | ~5 d / 14841 ms | — | OK |
| hd3.6 wp50 | 48 | 50 d / 4412 ms | ~5 d / 1746 ms | — | — | — | OK |
| hd3.7 machine | 10 | 15 d / 56 ms | ~7 d / 59 ms | 11 d / 111 ms | 15 d / 2551 ms | — | OK |
| hd3.7 wp50 | 48 | 50 d / 578 ms | ~7 d / 564 ms | — | — | — | OK |
| hd3.8 machine | 10 | 14 d / 51 ms | ~3 d / 105 ms | 13 d / 204 ms | 15 d / 4800 ms | — | FEWER-DIGITS(14<15) |
| hd3.8 wp50 | 48 | 50 d / 943 ms | ~3 d / 1211 ms | — | — | — | OK |
| hd3.9 machine | 10 | 14 d / 1244 ms | ~7 d / 544 ms | 11 d / 1695 ms | ~8 d / 100678 ms | — | OK |
| hd3.9 wp50 | 48 | 50 d / 11352 ms | ~7 d / 5587 ms | — | — | — | OK |
| hd3.10 machine | 10 | 14 d / 3754 ms | ~4 d / 1216 ms | 11 d / 5918 ms | — | — | OK |
| hd3.10 wp50 | 48 | 50 d / 36240 ms | ~4 d / 11999 ms | — | — | — | OK |

### 209-siam100h-d4-number-theory-roots-family

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd4.1 machine | 10 | 15 d / 0.2 ms | 15 d / 16 ms | 15 d / 1.3 ms | 15 d / 7 ms | — | OK |
| hd4.1 wp50 | 48 | 50 d / 1.4 ms | 50 d / 4.9 ms | 50 d / 1.8 ms | 50 d / 143 ms | — | OK |
| hd4.2 machine | 10 | 15 d / 0.07 ms | 15 d / 26 ms | 15 d / 7.6 ms | — | — | OK |
| hd4.2 wp50 | 48 | 50 d / 14 ms | 50 d / 73 ms | 50 d / 20 ms | — | — | OK |
| hd4.3 machine | 10 | 15 d / 0.25 ms | 15 d / 2646 ms | 14 d / 14 ms | — | — | OK |
| hd4.3 wp50 | 48 | 50 d / 19 ms | ~35 d / 3684 ms | 50 d / 72 ms | — | — | OK |
| hd4.4 machine | 10 | ~4 d / 10 ms | — | ~4 d / 73 ms | — | — | OK |
| hd4.4 wp50 | 48 | ~4 d / 252 ms | — | ~4 d / 112 ms | — | — | OK |
| hd4.5 machine | 10 | 15 d / 0.32 ms | 15 d / 3725 ms | 15 d / 17 ms | — | — | OK |
| hd4.5 wp50 | 48 | 50 d / 22 ms | ~34 d / 5019 ms | 50 d / 75 ms | — | — | OK |
| hd4.6 machine | 10 | 15 d / 0.14 ms | ~2 d / 18158 ms | ~1 d / 91 ms | — | — | OK |
| hd4.6 wp50 | 48 | 50 d / 2.1 ms | ~2 d / 33068 ms | ~1 d / 127 ms | — | — | OK |
| hd4.7 machine | 10 | 15 d / 1.3 ms | ~1 d / 33920 ms | ✗0 d / 7.2 ms | — | — | OK |
| hd4.7 wp50 | 48 | 50 d / 6.1 ms | — | ✗0 d / 7.8 ms | — | — | OK |
| hd4.8 machine | 10 | 15 d / 1.1 ms | 15 d / 89 ms | 15 d / 12 ms | — | — | OK |
| hd4.8 wp50 | 48 | 50 d / 4.6 ms | 50 d / 162 ms | 50 d / 25 ms | — | — | OK |
| hd4.9 machine | 10 | 14 d / 0.46 ms | ✗0 d / 37712 ms | 14 d / 350 ms | — | — | OK |
| hd4.9 wp50 | 48 | 50 d / 48 ms | — | 50 d / 982 ms | — | — | OK |
| hd4.10 machine | 10 | 11 d / 16049 ms | — | — | — | — | OK |
| hd4.10 wp50 | 48 | 50 d / 16894 ms | 50 d / 6602 ms | — | — | — | SLOWER(2.6x) |

### 210-siam100h-d5-matrix-polynomial-family

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd5.1 machine | 10 | ✗0 d / 279 ms | — | ✗0 d / 36 ms | — | — | OK(nobody) |
| hd5.1 wp50 | 48 | ~6 d / 3544 ms | — | — | — | — | OK(nobody) |
| hd5.2 machine | 10 | ✗0 d / 31 ms | ✗0 d / 235 ms | ✗0 d / 101 ms | — | — | OK(nobody) |
| hd5.2 wp50 | 48 | ✗0 d / 1860 ms | ✗0 d / 629 ms | — | — | — | OK(nobody) |
| hd5.3 machine | 10 | ~1 d / 1549 ms | 13 d / 5019 ms | 10 d / 50 ms | 13 d / 23 ms | — | FEWER-DIGITS(1<13) |
| hd5.3 wp50 | 48 | — | 50 d / 11986 ms | — | — | — | FEWER-DIGITS(0<50) |
| hd5.4 machine | 10 | ✗0 d / 30 ms | ✗0 d / 759 ms | ✗0 d / 134 ms | — | — | OK(nobody) |
| hd5.4 wp50 | 48 | ~1 d / 18672 ms | ✗0 d / 41069 ms | — | — | — | OK(nobody) |
| hd5.5 machine | 10 | ~3 d / 548 ms | ✗0 d / 2137 ms | ~1 d / 146 ms | — | — | OK |
| hd5.5 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hd5.6 machine | 10 | ~4 d / 579 ms | ✗0 d / 918 ms | ~3 d / 340 ms | — | — | OK |
| hd5.6 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hd5.7 machine | 10 | ~1 d / 4 ms | ✗0 d / 49 ms | ✗0 d / 43 ms | — | — | OK(nobody) |
| hd5.7 wp50 | 48 | ~2 d / 1076 ms | ✗0 d / 273 ms | — | — | — | OK(nobody) |
| hd5.8 machine | 10 | ~5 d / 12 ms | ✗0 d / 215 ms | ~7 d / 91 ms | — | — | FEWER-DIGITS(5<7) |
| hd5.8 wp50 | 48 | ~19 d / 2785 ms | ✗0 d / 8034 ms | — | — | — | OK(nobody) |
| hd5.9 machine | 10 | 15 d / 678 ms | ✗0 d / 265 ms | ✗0 d / 36 ms | — | — | OK |
| hd5.9 wp50 | 48 | 50 d / 6491 ms | ✗0 d / 13616 ms | — | — | — | OK |
| hd5.10 machine | 10 | ✗0 d / 5553 ms | ✗0 d / 13303 ms | ✗0 d / 679 ms | — | — | OK(nobody) |
| hd5.10 wp50 | 48 | — | — | — | — | — | OK(nobody) |

### 211-siam100h-d6-oscillatory-towers

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd6.1 machine | 10 | 13 d / 0.14 ms | 13 d / 5.9 ms | 13 d / 0.2 ms | 13 d / 7.5 ms | 13 d / 11 ms | OK |
| hd6.1 wp50 | 48 | 50 d / 98 ms | 50 d / 107 ms | 50 d / 54 ms | ~7 d / 276 ms | 50 d / 29440 ms | SLOWER(1.8x) |
| hd6.2 machine | 10 | 15 d / 0.063 ms | 12 d / 6 ms | 15 d / 0.092 ms | 15 d / 3.7 ms | 14 d / 3.9 ms | OK |
| hd6.2 wp50 | 48 | 50 d / 27 ms | 50 d / 62 ms | 50 d / 48 ms | ~7 d / 134 ms | 50 d / 26325 ms | OK |
| hd6.3 machine | 10 | 15 d / 1.2 ms | 13 d / 7.8 ms | 13 d / 2.7 ms | 13 d / 40 ms | 14 d / 72 ms | OK |
| hd6.3 wp50 | 48 | 50 d / 217 ms | 50 d / 926 ms | ✗0 d / 100 ms | — | ~24 d / 37215 ms | OK |
| hd6.4 machine | 10 | 15 d / 0.045 ms | 15 d / 7.2 ms | 15 d / 0.067 ms | 15 d / 3.8 ms | 15 d / 2.5 ms | OK |
| hd6.4 wp50 | 48 | 50 d / 77 ms | 50 d / 95 ms | 50 d / 13 ms | ~7 d / 17 ms | 50 d / 2941 ms | SLOWER(6.0x) |
| hd6.5 machine | 10 | 15 d / 0.061 ms | 12 d / 4.5 ms | 15 d / 0.097 ms | 15 d / 4.2 ms | 13 d / 4.4 ms | OK |
| hd6.5 wp50 | 48 | 50 d / 58 ms | 50 d / 71 ms | 50 d / 93 ms | ~6 d / 54 ms | 50 d / 19034 ms | OK |
| hd6.6 machine | 10 | 15 d / 0.14 ms | 13 d / 6.9 ms | 15 d / 0.14 ms | 15 d / 6.1 ms | 14 d / 4.1 ms | OK |
| hd6.6 wp50 | 48 | 50 d / 110 ms | 50 d / 66 ms | 50 d / 37 ms | ~7 d / 131 ms | 50 d / 25383 ms | SLOWER(3.0x) |
| hd6.7 machine | 10 | 13 d / 0.15 ms | 12 d / 16 ms | 13 d / 0.29 ms | 13 d / 8.4 ms | 14 d / 6.9 ms | FEWER-DIGITS(13<14) |
| hd6.7 wp50 | 48 | 50 d / 0.52 ms | 50 d / 139 ms | 50 d / 43 ms | ~8 d / 373 ms | 50 d / 35704 ms | OK |
| hd6.8 machine | 10 | 15 d / 0.33 ms | 11 d / 24 ms | 15 d / 0.48 ms | 15 d / 14 ms | 15 d / 74 ms | OK |
| hd6.8 wp50 | 48 | 50 d / 202 ms | 50 d / 311 ms | ~18 d / 120 ms | ~7 d / 473 ms | ~42 d / 47846 ms | OK |
| hd6.9 machine | 10 | 15 d / 0.075 ms | 14 d / 9.9 ms | 14 d / 0.18 ms | 14 d / 6.5 ms | 14 d / 2.8 ms | OK |
| hd6.9 wp50 | 48 | 50 d / 82 ms | 50 d / 65 ms | 50 d / 24 ms | ~7 d / 26 ms | 50 d / 3526 ms | SLOWER(3.4x) |
| hd6.10 machine | 10 | 15 d / 0.14 ms | 14 d / 13 ms | 15 d / 0.18 ms | 15 d / 5.2 ms | 13 d / 505 ms | OK |
| hd6.10 wp50 | 48 | 50 d / 107 ms | 50 d / 143 ms | 50 d / 69 ms | — | ~31 d / 20100 ms | SLOWER(1.5x) |

### 212-siam100h-d7-region-areas

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd7.1 machine | 10 | 15 d / 0.13 ms | 14 d / 6.6 ms | 15 d / 0.17 ms | 15 d / 4.3 ms | ~8 d / 1156 ms | OK |
| hd7.1 wp50 | 48 | 50 d / 8.9 ms | 50 d / 77 ms | 50 d / 36 ms | — | ~7 d / 6981 ms | OK |
| hd7.2 machine | 10 | 15 d / 0.069 ms | 15 d / 8 ms | 15 d / 0.22 ms | 15 d / 9 ms | — | OK |
| hd7.2 wp50 | 48 | 50 d / 1.4 ms | 50 d / 24 ms | 50 d / 6.6 ms | ~6 d / 13 ms | — | OK |
| hd7.3 machine | 10 | 15 d / 6.9 ms | 14 d / 25 ms | 14 d / 4.2 ms | 15 d / 122 ms | — | SLOWER(1.6x) |
| hd7.3 wp50 | 48 | 50 d / 62 ms | 50 d / 764 ms | 50 d / 88 ms | ~4 d / 3529 ms | — | OK |
| hd7.4 machine | 10 | 15 d / 0.5 ms | 13 d / 16 ms | 15 d / 3.2 ms | 15 d / 51 ms | — | OK |
| hd7.4 wp50 | 48 | 50 d / 65 ms | 50 d / 466 ms | 50 d / 87 ms | ~4 d / 1494 ms | — | OK |
| hd7.5 machine | 10 | 15 d / 6.8 ms | 14 d / 5048 ms | 15 d / 7.2 ms | 14 d / 148 ms | — | OK |
| hd7.5 wp50 | 48 | 50 d / 95 ms | 50 d / 6083 ms | 50 d / 133 ms | — | — | OK |
| hd7.6 machine | 10 | 15 d / 7.3 ms | 14 d / 5051 ms | 14 d / 5.7 ms | 15 d / 648 ms | — | OK |
| hd7.6 wp50 | 48 | 50 d / 94 ms | 50 d / 6081 ms | 50 d / 108 ms | ~4 d / 5426 ms | — | OK |
| hd7.7 machine | 10 | 15 d / 0.074 ms | 13 d / 15 ms | 15 d / 0.17 ms | 15 d / 7.4 ms | 10 d / 822 ms | OK |
| hd7.7 wp50 | 48 | 50 d / 23 ms | 50 d / 199 ms | 50 d / 33 ms | ~4 d / 135 ms | ~9 d / 46737 ms | OK |
| hd7.8 machine | 10 | 15 d / 0.053 ms | 12 d / 5.5 ms | 15 d / 0.063 ms | 15 d / 4.2 ms | 14 d / 3.5 ms | OK |
| hd7.8 wp50 | 48 | 50 d / 13 ms | 50 d / 58 ms | 50 d / 15 ms | ~6 d / 28 ms | 50 d / 26936 ms | OK |
| hd7.9 machine | 10 | 15 d / 0.1 ms | 15 d / 8.5 ms | 15 d / 0.38 ms | 15 d / 7.8 ms | — | OK |
| hd7.9 wp50 | 48 | 50 d / 1.2 ms | 50 d / 23 ms | 50 d / 6.5 ms | ~9 d / 17 ms | — | OK |
| hd7.10 machine | 10 | 15 d / 7 ms | 13 d / 9976 ms | 14 d / 7.1 ms | 15 d / 166 ms | — | OK |
| hd7.10 wp50 | 48 | 50 d / 149 ms | 50 d / 10724 ms | 50 d / 176 ms | ~5 d / 5882 ms | — | OK |

### 213-siam100h-d8-sobolev-family

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd8.1 machine | 10 | 15 d / 0.024 ms | 15 d / 5.2 ms | 15 d / 0.053 ms | 15 d / 4.4 ms | — | OK |
| hd8.1 wp50 | 48 | 50 d / 2.3 ms | 50 d / 20 ms | 50 d / 27 ms | — | — | OK |
| hd8.2 machine | 10 | 15 d / 0.04 ms | 10 d / 5.6 ms | 15 d / 0.092 ms | 15 d / 2.7 ms | — | OK |
| hd8.2 wp50 | 48 | 50 d / 7.5 ms | 50 d / 30 ms | 50 d / 46 ms | — | — | OK |
| hd8.3 machine | 10 | 15 d / 0.064 ms | 14 d / 3.9 ms | 15 d / 0.28 ms | 15 d / 3.5 ms | — | OK |
| hd8.3 wp50 | 48 | 50 d / 12 ms | 50 d / 38 ms | 50 d / 177 ms | — | — | OK |
| hd8.4 machine | 10 | 15 d / 0.13 ms | 11 d / 6.3 ms | 15 d / 0.56 ms | 15 d / 4.8 ms | — | OK |
| hd8.4 wp50 | 48 | 50 d / 23 ms | 50 d / 96 ms | 50 d / 287 ms | — | — | OK |
| hd8.5 machine | 10 | 15 d / 30 ms | ~5 d / 32 ms | 15 d / 200 ms | 15 d / 1124 ms | — | OK |
| hd8.5 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hd8.6 machine | 10 | 15 d / 49 ms | ~8 d / 18 ms | 15 d / 532 ms | 15 d / 1983 ms | — | OK |
| hd8.6 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hd8.7 machine | 10 | — | ~8 d / 3167 ms | — | — | — | FEWER-DIGITS(0<8) |
| hd8.7 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hd8.8 machine | 10 | 15 d / 83 ms | ~7 d / 48 ms | 15 d / 19 ms | 15 d / 576 ms | — | SLOWER(4.3x) |
| hd8.8 wp50 | 48 | 50 d / 30789 ms | — | — | — | — | OK |
| hd8.9 machine | 10 | — | ~8 d / 1599 ms | — | — | — | FEWER-DIGITS(0<8) |
| hd8.9 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hd8.10 machine | 10 | 15 d / 36 ms | ~7 d / 190 ms | 15 d / 110 ms | 15 d / 1093 ms | — | OK |
| hd8.10 wp50 | 48 | 50 d / 53787 ms | — | — | — | — | OK |

### 214-siam100h-d9-lattice-walks

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd9.1 machine | 10 | 15 d / 0.065 ms | 15 d / 0.19 ms | 15 d / 0.098 ms | 15 d / 219 ms | 15 d / 61 ms | OK |
| hd9.1 wp50 | 48 | 50 d / 4 ms | 50 d / 7 ms | 50 d / 48 ms | 50 d / 13 ms | 50 d / 17 ms | OK |
| hd9.2 machine | 10 | 15 d / 7.6 ms | 14 d / 18 ms | 15 d / 12 ms | — | 14 d / 41367 ms | OK |
| hd9.2 wp50 | 48 | 50 d / 624 ms | 50 d / 1028 ms | 50 d / 9027 ms | 50 d / 2204 ms | 50 d / 1492 ms | OK |
| hd9.3 machine | 10 | 14 d / 1.5 ms | 14 d / 2.7 ms | 14 d / 1.5 ms | 14 d / 21275 ms | 15 d / 3225 ms | FEWER-DIGITS(14<15) |
| hd9.3 wp50 | 48 | 50 d / 118 ms | 50 d / 306 ms | 50 d / 1465 ms | 50 d / 511 ms | 50 d / 362 ms | OK |
| hd9.4 machine | 10 | 14 d / 5.6 ms | 14 d / 12 ms | 15 d / 5.7 ms | — | 15 d / 26876 ms | FEWER-DIGITS(14<15) |
| hd9.4 wp50 | 48 | 50 d / 527 ms | 50 d / 605 ms | 50 d / 6674 ms | 50 d / 1000 ms | 50 d / 810 ms | OK |
| hd9.5 machine | 10 | 15 d / 0.76 ms | 15 d / 1 ms | 15 d / 0.86 ms | 15 d / 3921 ms | 15 d / 1650 ms | OK |
| hd9.5 wp50 | 48 | 50 d / 34 ms | 50 d / 134 ms | 50 d / 424 ms | 50 d / 146 ms | 50 d / 166 ms | OK |
| hd9.6 machine | 10 | 14 d / 1.5 ms | 15 d / 2.3 ms | 14 d / 1.5 ms | 15 d / 18100 ms | 15 d / 4666 ms | FEWER-DIGITS(14<15) |
| hd9.6 wp50 | 48 | 50 d / 96 ms | 50 d / 168 ms | 50 d / 1280 ms | 50 d / 241 ms | 50 d / 237 ms | OK |
| hd9.7 machine | 10 | 15 d / 2.1 ms | 15 d / 4.2 ms | 15 d / 2.2 ms | 15 d / 21758 ms | 14 d / 6804 ms | OK |
| hd9.7 wp50 | 48 | 50 d / 102 ms | 50 d / 176 ms | 50 d / 1333 ms | 50 d / 240 ms | 50 d / 129 ms | OK |
| hd9.8 machine | 10 | 14 d / 7.7 ms | 14 d / 18 ms | 14 d / 11 ms | — | 14 d / 40115 ms | OK |
| hd9.8 wp50 | 48 | 50 d / 631 ms | 50 d / 1152 ms | 50 d / 8860 ms | 50 d / 2269 ms | 50 d / 1526 ms | OK |
| hd9.9 machine | 10 | 14 d / 7.8 ms | 14 d / 17 ms | 14 d / 11 ms | — | 15 d / 37799 ms | FEWER-DIGITS(14<15) |
| hd9.9 wp50 | 48 | 50 d / 630 ms | 50 d / 1004 ms | 50 d / 8693 ms | 50 d / 2271 ms | 50 d / 1323 ms | OK |
| hd9.10 machine | 10 | 14 d / 32 ms | 14 d / 67 ms | 14 d / 36 ms | — | 14 d / 185877 ms | OK |
| hd9.10 wp50 | 48 | 50 d / 3038 ms | 50 d / 3688 ms | 50 d / 36371 ms | 50 d / 10177 ms | 50 d / 4279 ms | OK |

### 215-siam100h-d10-random-recurrences

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd10.1 machine | 2 | 3 d / 7.8 ms | 5 d / 4.6 ms | 4 d / 169 ms | 5 d / 18994 ms | 3 d / 165 ms | FEWER-DIGITS(3<5) SLOWER(1.7x) |
| hd10.1 wp50 | 1 | 3 d / 5.4 ms | 2 d / 3.9 ms | 3 d / 2.1 ms | 3 d / 149 ms | 3 d / 25 ms | SLOWER(2.5x) |
| hd10.2 machine | 2 | 4 d / 24 ms | 5 d / 127 ms | 4 d / 222 ms | 4 d / 33205 ms | 3 d / 417 ms | FEWER-DIGITS(4<5) |
| hd10.2 wp50 | 1 | 3 d / 27 ms | 2 d / 37 ms | 3 d / 14 ms | 3 d / 5509 ms | 3 d / 28 ms | SLOWER(2.0x) |
| hd10.3 machine | 2 | 4 d / 24 ms | 5 d / 146 ms | 4 d / 216 ms | 5 d / 18903 ms | 4 d / 354 ms | FEWER-DIGITS(4<5) |
| hd10.3 wp50 | 1 | 4 d / 28 ms | 2 d / 50 ms | 3 d / 16 ms | 3 d / 13254 ms | 3 d / 12 ms | SLOWER(2.4x) |
| hd10.4 machine | 2 | 4 d / 24 ms | 4 d / 155 ms | 4 d / 446 ms | 4 d / 19345 ms | 3 d / 2675 ms | OK |
| hd10.4 wp50 | 1 | 3 d / 258 ms | 2 d / 180 ms | 3 d / 3796 ms | 3 d / 19852 ms | 2 d / 401 ms | OK |
| hd10.5 machine | 2 | 4 d / 26 ms | 4 d / 163 ms | 4 d / 588 ms | 4 d / 20405 ms | 4 d / 330 ms | OK |
| hd10.5 wp50 | 1 | 3 d / 27 ms | 2 d / 51 ms | 3 d / 31 ms | 2 d / 6505 ms | 2 d / 9.8 ms | SLOWER(2.8x) |
| hd10.6 machine | 2 | 3 d / 20 ms | 5 d / 179 ms | 5 d / 652 ms | 5 d / 16855 ms | 3 d / 111 ms | FEWER-DIGITS(3<5) |
| hd10.6 wp50 | 1 | 3 d / 21 ms | 2 d / 31 ms | 4 d / 12 ms | 3 d / 1971 ms | 2 d / 5.2 ms | FEWER-DIGITS(3<4) SLOWER(4.0x) |
| hd10.7 machine | 2 | 4 d / 24 ms | 4 d / 172 ms | 4 d / 560 ms | 4 d / 17248 ms | 4 d / 338 ms | OK |
| hd10.7 wp50 | 1 | 3 d / 27 ms | 3 d / 34 ms | 3 d / 20 ms | 2 d / 226 ms | 2 d / 10 ms | SLOWER(2.7x) |
| hd10.8 machine | 2 | 4 d / 28 ms | 4 d / 168 ms | 4 d / 581 ms | 5 d / 14724 ms | 4 d / 287 ms | FEWER-DIGITS(4<5) |
| hd10.8 wp50 | 1 | 3 d / 28 ms | 2 d / 45 ms | 5 d / 24 ms | 3 d / 6776 ms | 3 d / 9.6 ms | FEWER-DIGITS(3<5) SLOWER(2.9x) |
| hd10.9 machine | 2 | 4 d / 21 ms | 4 d / 195 ms | 3 d / 608 ms | 3 d / 19044 ms | 3 d / 118 ms | OK |
| hd10.9 wp50 | 1 | 2 d / 22 ms | 2 d / 47 ms | 2 d / 14 ms | 2 d / 3200 ms | 2 d / 7.7 ms | SLOWER(2.8x) |
| hd10.10 machine | 2 | 5 d / 28 ms | 3 d / 150 ms | 3 d / 562 ms | 5 d / 15048 ms | 3 d / 339 ms | OK |
| hd10.10 wp50 | 1 | 3 d / 26 ms | 2 d / 35 ms | 3 d / 28 ms | 2 d / 2642 ms | 2 d / 12 ms | SLOWER(2.2x) |

### 216-siam100h-d11-shooting-targets

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd11.1 machine | 10 | 15 d / 1.8 ms | 13 d / 45 ms | 15 d / 4.9 ms | — | — | OK |
| hd11.1 wp50 | 48 | 50 d / 14 ms | 50 d / 111 ms | 50 d / 804 ms | — | — | OK |
| hd11.2 machine | 10 | 13 d / 2.9 ms | 12 d / 23 ms | 13 d / 8.2 ms | — | — | OK |
| hd11.2 wp50 | 48 | 50 d / 17 ms | ~2 d / 10974 ms | 50 d / 2389 ms | — | — | OK |
| hd11.3 machine | 10 | 13 d / 2.4 ms | 13 d / 13 ms | 14 d / 8.7 ms | — | — | FEWER-DIGITS(13<14) |
| hd11.3 wp50 | 48 | 50 d / 13 ms | ~3 d / 9938 ms | 50 d / 2003 ms | — | — | OK |
| hd11.4 machine | 10 | 14 d / 1.9 ms | 13 d / 15 ms | 14 d / 4.1 ms | — | — | OK |
| hd11.4 wp50 | 48 | 50 d / 16 ms | ~2 d / 9582 ms | 50 d / 1991 ms | — | — | OK |
| hd11.5 machine | 10 | 12 d / 5.8 ms | 11 d / 41 ms | 12 d / 24 ms | — | — | OK |
| hd11.5 wp50 | 48 | 50 d / 692 ms | — | 50 d / 19746 ms | — | — | OK |
| hd11.6 machine | 10 | 13 d / 4.7 ms | 11 d / 29 ms | 12 d / 31 ms | — | — | OK |
| hd11.6 wp50 | 48 | 50 d / 531 ms | — | 50 d / 14050 ms | — | — | OK |
| hd11.7 machine | 10 | 13 d / 13 ms | 11 d / 76 ms | 13 d / 98 ms | — | — | OK |
| hd11.7 wp50 | 48 | 50 d / 3228 ms | — | 50 d / 74348 ms | — | — | OK |
| hd11.8 machine | 10 | 13 d / 16 ms | 12 d / 62 ms | 13 d / 140 ms | — | — | OK |
| hd11.8 wp50 | 48 | 50 d / 4345 ms | — | 50 d / 86846 ms | — | — | OK |
| hd11.9 machine | 10 | 13 d / 3 ms | 12 d / 12 ms | 13 d / 12 ms | — | — | OK |
| hd11.9 wp50 | 48 | 50 d / 15 ms | ~2 d / 5761 ms | 50 d / 4229 ms | — | — | OK |
| hd11.10 machine | 10 | 12 d / 93 ms | 11 d / 193 ms | 12 d / 1204 ms | — | — | OK |
| hd11.10 wp50 | 48 | 50 d / 30586 ms | — | — | — | — | OK |

### 217-siam100h-d12-event-parameter

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd12.1 machine | 10 | 12 d / 3.3 ms | 10 d / 18 ms | 12 d / 21 ms | 11 d / 224 ms | — | OK |
| hd12.1 wp50 | 48 | 50 d / 66 ms | ~46 d / 44231 ms | 50 d / 9356 ms | — | — | OK |
| hd12.2 machine | 10 | 12 d / 1.5 ms | 11 d / 13 ms | 12 d / 13 ms | 12 d / 180 ms | — | OK |
| hd12.2 wp50 | 48 | 50 d / 67 ms | ~46 d / 67953 ms | 50 d / 31645 ms | — | — | OK |
| hd12.3 machine | 10 | 12 d / 1.8 ms | 10 d / 12 ms | 12 d / 16 ms | 10 d / 112 ms | — | OK |
| hd12.3 wp50 | 48 | 50 d / 55 ms | — | 50 d / 4946 ms | — | — | OK |
| hd12.4 machine | 10 | 11 d / 1.5 ms | 10 d / 9.4 ms | 11 d / 13 ms | 13 d / 119 ms | — | FEWER-DIGITS(11<13) |
| hd12.4 wp50 | 48 | 50 d / 5 ms | — | 50 d / 2154 ms | — | — | OK |
| hd12.5 machine | 10 | 12 d / 0.89 ms | 10 d / 13 ms | 12 d / 25 ms | 12 d / 37 ms | — | OK |
| hd12.5 wp50 | 48 | 50 d / 19 ms | — | 50 d / 1942 ms | — | — | OK |
| hd12.6 machine | 10 | 11 d / 1.4 ms | 10 d / 9.4 ms | 12 d / 18 ms | 10 d / 87 ms | — | FEWER-DIGITS(11<12) |
| hd12.6 wp50 | 48 | 50 d / 136 ms | — | 50 d / 9871 ms | — | — | OK |
| hd12.7 machine | 10 | 13 d / 1.8 ms | 11 d / 12 ms | 12 d / 19 ms | 10 d / 83 ms | — | OK |
| hd12.7 wp50 | 48 | 50 d / 12 ms | — | 50 d / 1443 ms | — | — | OK |
| hd12.8 machine | 10 | 12 d / 1.4 ms | 11 d / 9.5 ms | 12 d / 12 ms | 11 d / 111 ms | — | OK |
| hd12.8 wp50 | 48 | 50 d / 21 ms | 48 d / 44372 ms | 50 d / 11328 ms | — | — | OK |
| hd12.9 machine | 10 | 13 d / 4 ms | 11 d / 13 ms | 12 d / 31 ms | 11 d / 729 ms | — | OK |
| hd12.9 wp50 | 48 | 50 d / 844 ms | — | — | — | — | OK |
| hd12.10 machine | 10 | 13 d / 2 ms | 11 d / 8.9 ms | 12 d / 20 ms | 11 d / 144 ms | — | OK |
| hd12.10 wp50 | 48 | 50 d / 93 ms | — | 50 d / 8609 ms | — | — | OK |

### 218-siam100h-d13-initial-value-shooting

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd13.1 machine | 10 | ~8 d / 2.7 ms | ~7 d / 8.2 ms | ~8 d / 25 ms | 11 d / 437 ms | — | FEWER-DIGITS(8<11) |
| hd13.1 wp50 | 48 | 50 d / 74 ms | — | 50 d / 3740 ms | — | — | OK |
| hd13.2 machine | 10 | 11 d / 3.1 ms | 12 d / 4.5 ms | 11 d / 16 ms | 12 d / 292 ms | — | FEWER-DIGITS(11<12) |
| hd13.2 wp50 | 48 | 50 d / 6.2 ms | — | 50 d / 2497 ms | — | — | OK |
| hd13.3 machine | 10 | 12 d / 1.6 ms | 10 d / 6.3 ms | 12 d / 17 ms | 12 d / 135 ms | — | OK |
| hd13.3 wp50 | 48 | 50 d / 52 ms | — | 50 d / 2355 ms | — | — | OK |
| hd13.4 machine | 10 | 12 d / 1.1 ms | 10 d / 6.4 ms | 12 d / 7.4 ms | 13 d / 156 ms | — | FEWER-DIGITS(12<13) |
| hd13.4 wp50 | 48 | 50 d / 61 ms | ~46 d / 19030 ms | 50 d / 32429 ms | — | — | OK |
| hd13.5 machine | 10 | 12 d / 0.82 ms | 12 d / 6 ms | 12 d / 5.8 ms | 12 d / 81 ms | — | OK |
| hd13.5 wp50 | 48 | 50 d / 22 ms | ~46 d / 13656 ms | 50 d / 17386 ms | — | — | OK |
| hd13.6 machine | 10 | 12 d / 2.4 ms | 11 d / 6.3 ms | 12 d / 24 ms | 13 d / 508 ms | — | FEWER-DIGITS(12<13) |
| hd13.6 wp50 | 48 | 50 d / 71 ms | — | 50 d / 18076 ms | — | — | OK |
| hd13.7 machine | 10 | 12 d / 2 ms | 10 d / 16 ms | 12 d / 12 ms | 12 d / 289 ms | — | OK |
| hd13.7 wp50 | 48 | 50 d / 174 ms | ~47 d / 41777 ms | 50 d / 16830 ms | — | — | OK |
| hd13.8 machine | 10 | 12 d / 1.7 ms | 10 d / 9.4 ms | 12 d / 17 ms | 12 d / 73 ms | — | OK |
| hd13.8 wp50 | 48 | 50 d / 52 ms | ~45 d / 25484 ms | 50 d / 2106 ms | — | — | OK |
| hd13.9 machine | 10 | 12 d / 16 ms | 10 d / 58 ms | 11 d / 193 ms | — | — | OK |
| hd13.9 wp50 | 48 | 50 d / 2140 ms | — | — | — | — | OK |
| hd13.10 machine | 10 | 11 d / 6.2 ms | 10 d / 9.8 ms | 12 d / 62 ms | 11 d / 589 ms | — | FEWER-DIGITS(11<12) |
| hd13.10 wp50 | 48 | 50 d / 306 ms | — | 50 d / 27531 ms | — | — | OK |

### 219-siam100h-d14-elliptic-parameter

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd14.1 machine | 10 | 12 d / 22 ms | 11 d / 132 ms | 13 d / 351 ms | — | — | FEWER-DIGITS(12<13) |
| hd14.1 wp50 | 48 | 50 d / 73 ms | — | — | — | — | OK |
| hd14.2 machine | 10 | 11 d / 12 ms | 10 d / 94 ms | 12 d / 186 ms | — | — | FEWER-DIGITS(11<12) |
| hd14.2 wp50 | 48 | 50 d / 38 ms | — | 50 d / 30173 ms | — | — | OK |
| hd14.3 machine | 10 | 13 d / 29 ms | 11 d / 274 ms | 13 d / 235 ms | — | — | OK |
| hd14.3 wp50 | 48 | 50 d / 805 ms | — | — | — | — | OK |
| hd14.4 machine | 10 | ~8 d / 9926 ms | ~9 d / 14760 ms | — | — | — | FEWER-DIGITS(8<9) |
| hd14.4 wp50 | 48 | ~8 d / 12093 ms | ~9 d / 74723 ms | — | — | — | FEWER-DIGITS(8<9) |
| hd14.5 machine | 10 | ~7 d / 2580 ms | ~7 d / 4153 ms | — | — | — | OK |
| hd14.5 wp50 | 48 | ~7 d / 34731 ms | — | — | — | — | OK(nobody) |
| hd14.6 machine | 10 | ~9 d / 862 ms | ~3 d / 3173 ms | — | — | — | OK |
| hd14.6 wp50 | 48 | ~8 d / 8254 ms | ~3 d / 22916 ms | — | — | — | OK |
| hd14.7 machine | 10 | ~7 d / 370 ms | ~7 d / 8002 ms | — | — | — | OK |
| hd14.7 wp50 | 48 | ~7 d / 3579 ms | — | — | — | — | OK(nobody) |
| hd14.8 machine | 10 | ~7 d / 2840 ms | ~7 d / 7978 ms | — | — | — | OK |
| hd14.8 wp50 | 48 | ~7 d / 36488 ms | — | — | — | — | OK(nobody) |
| hd14.9 machine | 10 | ~6 d / 2397 ms | ~6 d / 5270 ms | — | — | — | OK |
| hd14.9 wp50 | 48 | ~6 d / 27015 ms | — | — | — | — | OK(nobody) |
| hd14.10 machine | 10 | ~8 d / 2360 ms | ~3 d / 5401 ms | — | — | — | OK |
| hd14.10 wp50 | 48 | ~8 d / 2794 ms | ~3 d / 40319 ms | — | — | — | OK |

### 220-siam100h-d15-blowup-times

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd15.1 machine | 10 | 12 d / 6.2 ms | 10 d / 23 ms | 13 d / 42 ms | — | — | FEWER-DIGITS(12<13) |
| hd15.1 wp50 | 48 | ~32 d / 80 ms | ✗0 d / 48 ms | — | — | — | OK(nobody) |
| hd15.2 machine | 10 | 12 d / 3.1 ms | 10 d / 3.2 ms | 13 d / 40 ms | — | — | FEWER-DIGITS(12<13) |
| hd15.2 wp50 | 48 | ~32 d / 85 ms | ✗0 d / 4.1 ms | — | — | — | OK(nobody) |
| hd15.3 machine | 10 | 11 d / 6.6 ms | 10 d / 6.3 ms | 13 d / 74 ms | — | — | FEWER-DIGITS(11<13) |
| hd15.3 wp50 | 48 | ~33 d / 108 ms | ~2 d / 672 ms | — | — | — | OK |
| hd15.4 machine | 10 | 12 d / 1 ms | 10 d / 2 ms | 13 d / 46 ms | — | — | FEWER-DIGITS(12<13) |
| hd15.4 wp50 | 48 | ~32 d / 80 ms | ✗0 d / 5 ms | — | — | — | OK(nobody) |
| hd15.5 machine | 10 | 11 d / 3.3 ms | 10 d / 3.2 ms | 12 d / 45 ms | — | — | FEWER-DIGITS(11<12) |
| hd15.5 wp50 | 48 | ~32 d / 208 ms | ~3 d / 638 ms | — | — | — | OK |
| hd15.6 machine | 10 | 12 d / 3.6 ms | 10 d / 3.1 ms | 13 d / 41 ms | — | — | FEWER-DIGITS(12<13) |
| hd15.6 wp50 | 48 | ~32 d / 95 ms | ✗0 d / 22 ms | — | — | — | OK(nobody) |
| hd15.7 machine | 10 | ~8 d / 151 ms | ~9 d / 128 ms | ~7 d / 843 ms | — | — | FEWER-DIGITS(8<9) |
| hd15.7 wp50 | 48 | ~33 d / 2315 ms | ~16 d / 14879 ms | — | — | — | OK |
| hd15.8 machine | 10 | ~7 d / 3804 ms | ~8 d / 163 ms | ~8 d / 1588 ms | — | — | FEWER-DIGITS(7<8) |
| hd15.8 wp50 | 48 | ~7 d / 5388 ms | ~9 d / 18241 ms | — | — | — | FEWER-DIGITS(7<9) |
| hd15.9 machine | 10 | ~5 d / 5.8 ms | ~5 d / 65 ms | ~4 d / 114 ms | — | — | OK |
| hd15.9 wp50 | 48 | ~5 d / 211 ms | ~5 d / 12283 ms | — | — | — | OK |
| hd15.10 machine | 10 | 12 d / 2.9 ms | 10 d / 3.2 ms | 12 d / 42 ms | — | — | OK |
| hd15.10 wp50 | 48 | ~32 d / 108 ms | ~1 d / 603 ms | — | — | — | OK |

### 221-siam100h-d16-heat-crossings

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd16.1 machine | 10 | ~7 d / 24 ms | ~3 d / 423 ms | — | — | — | OK |
| hd16.1 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hd16.2 machine | 10 | ~6 d / 224 ms | ~3 d / 863 ms | — | — | — | OK |
| hd16.2 wp50 | 48 | ~2 d / 43 ms | — | — | — | — | OK(nobody) |
| hd16.3 machine | 10 | ~5 d / 4.2 ms | ✗0 d / 163 ms | — | — | — | OK(nobody) |
| hd16.3 wp50 | 48 | ~5 d / 43 ms | ✗0 d / 2622 ms | — | — | — | OK(nobody) |
| hd16.4 machine | 10 | ~6 d / 3.4 ms | ~3 d / 5.5 ms | — | — | — | OK |
| hd16.4 wp50 | 48 | ~6 d / 38 ms | ~6 d / 2540 ms | — | — | — | OK |
| hd16.5 machine | 10 | ~6 d / 980 ms | ~3 d / 1475 ms | — | — | — | OK |
| hd16.5 wp50 | 48 | ✗0 d / 43 ms | — | — | — | — | OK(nobody) |
| hd16.6 machine | 10 | ~3 d / 7.8 ms | ✗0 d / 5.9 ms | — | — | — | OK(nobody) |
| hd16.6 wp50 | 48 | ~3 d / 43 ms | ~1 d / 409 ms | — | — | — | OK |
| hd16.7 machine | 10 | ~6 d / 232 ms | ~3 d / 2721 ms | — | — | — | OK |
| hd16.7 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hd16.8 machine | 10 | ~4 d / 987 ms | ~3 d / 3277 ms | — | — | — | OK |
| hd16.8 wp50 | 48 | — | — | — | — | — | OK(nobody) |
| hd16.9 machine | 10 | ~6 d / 18 ms | ✗0 d / 13 ms | — | — | — | OK(nobody) |
| hd16.9 wp50 | 48 | ~5 d / 45 ms | ✗0 d / 2574 ms | — | — | — | OK(nobody) |
| hd16.10 machine | 10 | ~4 d / 152 ms | ✗0 d / 54 ms | — | — | — | OK(nobody) |
| hd16.10 wp50 | 48 | ~4 d / 104 ms | ✗0 d / 38799 ms | — | — | — | OK(nobody) |

### 222-siam100h-d17-nodal-lines

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd17.1 machine | 10 | ~4 d / 410 ms | ✗0 d / 1585 ms | — | — | — | OK(nobody) |
| hd17.1 wp50 | 48 | ~4 d / 408 ms | ✗0 d / 182 ms | — | — | — | OK(nobody) |
| hd17.2 machine | 10 | ~7 d / 288 ms | ✗0 d / 119 ms | — | — | — | OK(nobody) |
| hd17.2 wp50 | 48 | ~7 d / 299 ms | ✗0 d / 121 ms | — | — | — | OK(nobody) |
| hd17.3 machine | 10 | ~1 d / 381 ms | ~1 d / 762 ms | — | — | — | OK |
| hd17.3 wp50 | 48 | ~1 d / 379 ms | ~1 d / 669 ms | — | — | — | OK |
| hd17.4 machine | 10 | ~6 d / 433 ms | ~1 d / 559 ms | — | — | — | OK |
| hd17.4 wp50 | 48 | ~6 d / 436 ms | ~1 d / 480 ms | — | — | — | OK |
| hd17.5 machine | 10 | ~4 d / 340 ms | ✗0 d / 392 ms | — | — | — | OK(nobody) |
| hd17.5 wp50 | 48 | ~4 d / 337 ms | ✗0 d / 391 ms | — | — | — | OK(nobody) |
| hd17.6 machine | 10 | ~3 d / 310 ms | ✗0 d / 360 ms | — | — | — | OK(nobody) |
| hd17.6 wp50 | 48 | ~3 d / 323 ms | ✗0 d / 368 ms | — | — | — | OK(nobody) |
| hd17.7 machine | 10 | ~4 d / 384 ms | ✗0 d / 402 ms | — | — | — | OK(nobody) |
| hd17.7 wp50 | 48 | ~4 d / 381 ms | ✗0 d / 404 ms | — | — | — | OK(nobody) |
| hd17.8 machine | 10 | ~2 d / 516 ms | ✗0 d / 171 ms | — | — | — | OK(nobody) |
| hd17.8 wp50 | 48 | ~2 d / 522 ms | ✗0 d / 168 ms | — | — | — | OK(nobody) |
| hd17.9 machine | 10 | ~6 d / 607 ms | ~1 d / 487 ms | — | — | — | OK |
| hd17.9 wp50 | 48 | ~6 d / 601 ms | ~1 d / 520 ms | — | — | — | OK |
| hd17.10 machine | 10 | ~2 d / 310 ms | ✗0 d / 296 ms | — | — | — | OK(nobody) |
| hd17.10 wp50 | 48 | ~2 d / 302 ms | ✗0 d / 288 ms | — | — | — | OK(nobody) |

### 223-siam100h-d18-dirichlet-eigenvalues

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd18.1 machine | 10 | 11 d / 372 ms | ~3 d / 416 ms | — | — | — | OK |
| hd18.1 wp50 | 48 | ~11 d / 363 ms | ~3 d / 137 ms | — | — | — | OK |
| hd18.2 machine | 10 | ~9 d / 381 ms | ✗0 d / 62 ms | — | — | — | OK(nobody) |
| hd18.2 wp50 | 48 | ~9 d / 386 ms | ✗0 d / 61 ms | — | — | — | OK(nobody) |
| hd18.3 machine | 10 | ~9 d / 281 ms | ~5 d / 64 ms | — | — | — | OK |
| hd18.3 wp50 | 48 | ~9 d / 284 ms | ~5 d / 65 ms | — | — | — | OK |
| hd18.4 machine | 10 | ~9 d / 286 ms | ~4 d / 44 ms | — | — | — | OK |
| hd18.4 wp50 | 48 | ~9 d / 300 ms | ~4 d / 44 ms | — | — | — | OK |
| hd18.5 machine | 10 | ~9 d / 292 ms | ~4 d / 44 ms | — | — | — | OK |
| hd18.5 wp50 | 48 | ~9 d / 290 ms | ~4 d / 45 ms | — | — | — | OK |
| hd18.6 machine | 10 | ~9 d / 318 ms | ~5 d / 94 ms | — | — | — | OK |
| hd18.6 wp50 | 48 | ~9 d / 327 ms | ~5 d / 92 ms | — | — | — | OK |
| hd18.7 machine | 10 | ~9 d / 296 ms | ~4 d / 89 ms | — | — | — | OK |
| hd18.7 wp50 | 48 | ~9 d / 298 ms | ~4 d / 83 ms | — | — | — | OK |
| hd18.8 machine | 10 | ~9 d / 294 ms | ~5 d / 76 ms | — | — | — | OK |
| hd18.8 wp50 | 48 | ~9 d / 298 ms | ~5 d / 82 ms | — | — | — | OK |
| hd18.9 machine | 10 | ~7 d / 294 ms | ~4 d / 90 ms | — | — | — | OK |
| hd18.9 wp50 | 48 | ~7 d / 290 ms | ~4 d / 91 ms | — | — | — | OK |
| hd18.10 machine | 10 | ~8 d / 326 ms | ~3 d / 92 ms | — | — | — | OK |
| hd18.10 wp50 | 48 | ~8 d / 326 ms | ~3 d / 87 ms | — | — | — | OK |

### 224-siam100h-d19-design-variants

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd19.1 machine | 10 | 15 d / 148 ms | ~4 d / 2896 ms | 13 d / 59 ms | 13 d / 3649 ms | — | SLOWER(2.5x) |
| hd19.1 wp50 | 48 | 50 d / 482 ms | ~9 d / 8361 ms | — | — | — | OK |
| hd19.2 machine | 10 | 15 d / 175 ms | ~3 d / 2616 ms | 11 d / 79 ms | 14 d / 2188 ms | — | SLOWER(2.2x) |
| hd19.2 wp50 | 48 | 50 d / 511 ms | ~10 d / 7680 ms | — | — | — | OK |
| hd19.3 machine | 10 | 14 d / 1907 ms | ~2 d / 32813 ms | ~8 d / 2166 ms | 13 d / 43931 ms | — | OK |
| hd19.3 wp50 | 48 | 50 d / 4347 ms | ~9 d / 97739 ms | — | — | — | OK |
| hd19.4 machine | 10 | 14 d / 798 ms | ~2 d / 10506 ms | ~4 d / 905 ms | 14 d / 3558 ms | — | OK |
| hd19.4 wp50 | 48 | 50 d / 1607 ms | ~9 d / 31994 ms | — | — | — | OK |
| hd19.5 machine | 10 | 13 d / 32 ms | ~5 d / 231 ms | 12 d / 39 ms | 15 d / 208 ms | — | FEWER-DIGITS(13<15) |
| hd19.5 wp50 | 48 | 50 d / 92 ms | ~16 d / 582 ms | — | — | — | OK |
| hd19.6 machine | 10 | 15 d / 113 ms | ~4 d / 1701 ms | 13 d / 78 ms | 14 d / 695 ms | — | OK |
| hd19.6 wp50 | 48 | 50 d / 258 ms | ~8 d / 4289 ms | — | — | — | OK |
| hd19.7 machine | 10 | 13 d / 434 ms | ~1 d / 2625 ms | ~4 d / 347 ms | 14 d / 1862 ms | — | FEWER-DIGITS(13<14) |
| hd19.7 wp50 | 48 | 50 d / 685 ms | ~8 d / 7988 ms | — | — | — | OK |
| hd19.8 machine | 10 | 12 d / 1982 ms | ~2 d / 25094 ms | ~4 d / 2142 ms | 12 d / 25832 ms | — | OK |
| hd19.8 wp50 | 48 | 50 d / 3831 ms | ~9 d / 70071 ms | — | — | — | OK |
| hd19.9 machine | 10 | 15 d / 53 ms | ~3 d / 2445 ms | 15 d / 25 ms | 15 d / 475 ms | — | SLOWER(2.1x) |
| hd19.9 wp50 | 48 | 50 d / 55 ms | ~8 d / 9004 ms | — | — | — | OK |
| hd19.10 machine | 10 | 13 d / 2310 ms | ~1 d / 33313 ms | ~1 d / 168 ms | 12 d / 74981 ms | — | OK |
| hd19.10 wp50 | 48 | 50 d / 5925 ms | ~11 d / 98672 ms | — | — | — | OK |

### 225-siam100h-d20-half-iterates

| case | need | Mathilda | WL | Py | Maxima | FriCAS | verdict |
|---|---|---|---|---|---|---|---|
| hd20.1 machine | 10 | 15 d / 9.9 ms | 15 d / 30 ms | 15 d / 3.8 ms | 15 d / 20 ms | 15 d / 294 ms | SLOWER(2.6x) |
| hd20.1 wp50 | 48 | 50 d / 123 ms | 50 d / 136 ms | 50 d / 1301 ms | 50 d / 388 ms | 50 d / 3724 ms | OK |
| hd20.2 machine | 10 | 15 d / 3.3 ms | 15 d / 10 ms | 15 d / 1.7 ms | 15 d / 7.7 ms | 15 d / 180 ms | SLOWER(1.9x) |
| hd20.2 wp50 | 48 | 50 d / 74 ms | 50 d / 176 ms | 50 d / 1338 ms | 50 d / 289 ms | 50 d / 2450 ms | OK |
| hd20.3 machine | 10 | 15 d / 2.7 ms | 15 d / 5.3 ms | 15 d / 1.5 ms | 15 d / 8 ms | 15 d / 188 ms | SLOWER(1.8x) |
| hd20.3 wp50 | 48 | 50 d / 42 ms | 50 d / 133 ms | 50 d / 1232 ms | 50 d / 271 ms | 50 d / 2239 ms | OK |
| hd20.4 machine | 10 | 15 d / 2.6 ms | 15 d / 5.4 ms | 15 d / 1.5 ms | 15 d / 7.1 ms | 15 d / 182 ms | SLOWER(1.8x) |
| hd20.4 wp50 | 48 | 50 d / 41 ms | 50 d / 125 ms | 50 d / 1232 ms | 50 d / 272 ms | 50 d / 2152 ms | OK |
| hd20.5 machine | 10 | 15 d / 3.6 ms | 15 d / 12 ms | 15 d / 1.8 ms | 15 d / 11 ms | 15 d / 183 ms | SLOWER(2.0x) |
| hd20.5 wp50 | 48 | 50 d / 48 ms | 50 d / 1022 ms | 50 d / 1248 ms | 50 d / 310 ms | 50 d / 2402 ms | OK |
| hd20.6 machine | 10 | 15 d / 3.9 ms | 15 d / 8.4 ms | 15 d / 2.1 ms | 15 d / 9 ms | — | SLOWER(1.8x) |
| hd20.6 wp50 | 48 | 50 d / 68 ms | 50 d / 163 ms | 50 d / 2440 ms | 50 d / 456 ms | — | OK |
| hd20.7 machine | 10 | 15 d / 9.1 ms | 15 d / 35 ms | 15 d / 4 ms | 15 d / 11 ms | — | SLOWER(2.3x) |
| hd20.7 wp50 | 48 | 50 d / 114 ms | 50 d / 2168 ms | 50 d / 1341 ms | 50 d / 453 ms | — | OK |
| hd20.8 machine | 10 | 11 d / 11 ms | 10 d / 72 ms | 11 d / 2.9 ms | 11 d / 11 ms | — | SLOWER(3.6x) |
| hd20.8 wp50 | 48 | 49 d / 268 ms | 50 d / 429 ms | 50 d / 2533 ms | ~45 d / 684 ms | — | FEWER-DIGITS(49<50) |
| hd20.9 machine | 10 | 15 d / 5.3 ms | 15 d / 13 ms | 15 d / 6.2 ms | 15 d / 24 ms | — | OK |
| hd20.9 wp50 | 48 | 50 d / 216 ms | 50 d / 1862 ms | 50 d / 2864 ms | 50 d / 2421 ms | — | OK |
| hd20.10 machine | 10 | 15 d / 1.8 ms | 15 d / 2 ms | 15 d / 0.28 ms | 15 d / 5 ms | — | SLOWER(6.4x) |
| hd20.10 wp50 | 48 | 50 d / 2.1 ms | 50 d / 3.5 ms | 50 d / 13 ms | 50 d / 33 ms | — | OK |
