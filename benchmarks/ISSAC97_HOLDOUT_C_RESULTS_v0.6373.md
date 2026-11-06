# ISSAC '97 holdout C (156-165): classification against the true answers

| family | cases | baseline correct / wrong / none | Mathilda correct / wrong / none | both correct | Mathilda faster |
|---|---|---|---|---|---|
| 156 c1 ill-conditioned linear algebra | 29 | 29 / 0 / 0 | 29 / 0 / 0 | 29 | 22 |
| 157 c2 hard definite integrals | 20 | 11 / 8 / 1 | 18 / 0 / 2 | 10 | 8 |
| 158 c3 slowly convergent series | 20 | 17 / 3 / 0 | 20 / 0 / 0 | 17 | 15 |
| 159 c4 coefficients and residues | 30 | 29 / 1 / 0 | 26 / 2 / 2 | 26 | 18 |
| 160 c5 zeros of special functions | 20 | 19 / 1 / 0 | 16 / 0 / 4 | 15 | 13 |
| 161 c6 polynomial systems | 30 | 30 / 0 / 0 | 30 / 0 / 0 | 30 | 29 |
| 162 c7 cancellation limits | 30 | 28 / 1 / 1 | 27 / 1 / 2 | 26 | 23 |
| 163 c8 closed-form integrals | 30 | 28 / 0 / 2 | 26 / 0 / 4 | 26 | 25 |
| 164 c9 operator eigenvalues | 20 | 20 / 0 / 0 | 20 / 0 / 0 | 20 | 17 |
| 165 c10 ODE singularities | 20 | 20 / 0 / 0 | 16 / 0 / 4 | 16 | 15 |
| **total** | 249 | 231 / 14 / 4 | 228 / 3 / 18 | 215 | 185 |

| case | baseline | ms | Mathilda | ms | faster |
|---|---|---|---|---|---|
| c1.1 Det exact n=30 | correct | 9.49 | correct | 1.98 | Mathilda 4.8x |
| c1.1 SingularValueList machine n=8 | correct | 0.017 | correct | 0.017 | baseline 1x |
| c1.1 SingularValueList wp300 n=60 | correct | 175 | correct | 238 | baseline 1.4x |
| c1.2 Det exact n=40 | correct | 1940 | correct | 333 | Mathilda 5.8x |
| c1.2 Eigenvalues machine n=30 | correct | 0.034 | correct | 0.044 | baseline 1.3x |
| c1.2 Eigenvalues wp60 n=120 | correct | 1779 | correct | 149 | Mathilda 12x |
| c1.3 Inverse exact n=30 | correct | 8.31 | correct | 1.23 | Mathilda 6.7x |
| c1.3 SingularValueList machine n=12 | correct | 0.023 | correct | 0.021 | Mathilda 1.1x |
| c1.3 SingularValueList wp250 n=150 | correct | 2112 | correct | 1686 | Mathilda 1.3x |
| c1.4 LinearSolve exact n=60 | correct | 7.13 | correct | 4.29 | Mathilda 1.7x |
| c1.4 Eigenvalues machine n=20 | correct | 0.031 | correct | 0.029 | Mathilda 1.1x |
| c1.4 Eigenvalues wp60 n=60 | correct | 186 | correct | 7.2 | Mathilda 26x |
| c1.5 SingularValueList machine n=32 | correct | 0.073 | correct | 0.071 | Mathilda 1x |
| c1.5 SingularValueList wp60 n=48 | correct | 34.6 | correct | 22.5 | Mathilda 1.5x |
| c1.6 Det exact n=50 | correct | 97.8 | correct | 49 | Mathilda 2x |
| c1.6 Eigenvalues machine n=50 | correct | 0.082 | correct | 0.127 | baseline 1.5x |
| c1.6 Eigenvalues wp60 n=50 | correct | 166 | correct | 13.1 | Mathilda 13x |
| c1.7 Det exact n=14 | correct | 8.8 | correct | 0.332 | Mathilda 27x |
| c1.7 SingularValueList machine n=8 | correct | 0.016 | correct | 0.018 | baseline 1.1x |
| c1.7 SingularValueList wp60 n=20 | correct | 4.84 | correct | 3.45 | Mathilda 1.4x |
| c1.8 LeastSquares exact n=20 | correct | 17.4 | correct | 3.39 | Mathilda 5.1x |
| c1.8 LeastSquares machine n=20 | correct | 0.066 | correct | 0.049 | Mathilda 1.3x |
| c1.8 LeastSquares wp80 n=20 | correct | 6.13 | correct | 0.757 | Mathilda 8.1x |
| c1.9 Det exact n=30 | correct | 10.9 | correct | 3.6 | Mathilda 3x |
| c1.9 SingularValueList machine n=6 | correct | 0.019 | correct | 0.016 | Mathilda 1.2x |
| c1.9 SingularValueList wp400 n=60 | correct | 213 | correct | 290 | baseline 1.4x |
| c1.10 LinearSolve exact n=50 | correct | 9.09 | correct | 3.34 | Mathilda 2.7x |
| c1.10 LinearSolve machine n=60 | correct | 0.061 | correct | 0.105 | baseline 1.7x |
| c1.10 LinearSolve wp60 n=100 | correct | 17.7 | correct | 2.86 | Mathilda 6.2x |
| c2.1 NIntegrate machine | wrong | 82.5 | correct | 28 |  |
| c2.1 NIntegrate wp40 | wrong | 825 | none | - |  |
| c2.2 NIntegrate machine | correct | 15.7 | correct | 0.087 | Mathilda 1.8e+02x |
| c2.2 NIntegrate wp40 | correct | 145 | correct | 5.98 | Mathilda 24x |
| c2.4 NIntegrate machine | wrong | 11 | correct | 6961 |  |
| c2.4 NIntegrate wp40 | wrong | 102 | correct | 9616 |  |
| c2.5 NIntegrate machine | correct | 13.6 | correct | 186 | baseline 14x |
| c2.5 NIntegrate wp40 | correct | 414 | none | - |  |
| c2.6 NIntegrate machine | correct | 20.5 | correct | 0.123 | Mathilda 1.7e+02x |
| c2.6 NIntegrate wp40 | correct | 40.5 | correct | 5.58 | Mathilda 7.3x |
| c2.7 NIntegrate machine | correct | 24.9 | correct | 28.9 | baseline 1.2x |
| c2.7 NIntegrate wp40 | correct | 1916 | correct | 924 | Mathilda 2.1x |
| c2.9 NIntegrate machine | correct | 14.6 | correct | 1.6 | Mathilda 9.1x |
| c2.9 NIntegrate wp40 | none | - | correct | 437 |  |
| c2.10 NIntegrate machine | correct | 17.6 | correct | 0.34 | Mathilda 52x |
| c2.10 NIntegrate wp40 | wrong | 100 | correct | 6904 |  |
| c2.3 NIntegrate machine | wrong | 3.81 | correct | 0.116 |  |
| c2.3 NIntegrate wp40 | wrong | 18.1 | correct | 8.92 |  |
| c2.8 NIntegrate machine | correct | 4.19 | correct | 0.042 | Mathilda 1e+02x |
| c2.8 NIntegrate wp40 | wrong | 92.7 | correct | 41.4 |  |
| c3.1 Sum primes machine | correct | 163 | correct | 7121 | baseline 44x |
| c3.1 MoebiusMu log-zeta wp50 | correct | 8.34 | correct | 1.01 | Mathilda 8.2x |
| c3.2 Sum machine | correct | 3543 | correct | 186 | Mathilda 19x |
| c3.2 Sum Infinity wp30 | correct | 31.9 | correct | 0.124 | Mathilda 2.6e+02x |
| c3.3 NSum machine | correct | 135 | correct | 0.201 | Mathilda 6.7e+02x |
| c3.3 NSum wp40 | wrong | 137 | correct | 1.63 |  |
| c3.4 NSum machine | wrong | 1.2 | correct | 4.74 |  |
| c3.4 NSum wp40 | correct | 104 | correct | 252 | baseline 2.4x |
| c3.5 NSum machine | correct | 11.7 | correct | 0.771 | Mathilda 15x |
| c3.5 NSum wp40 | correct | 50.3 | correct | 24.7 | Mathilda 2x |
| c3.6 NSum machine | correct | 36.9 | correct | 0.094 | Mathilda 3.9e+02x |
| c3.6 NSum wp40 | wrong | 37.3 | correct | 0.543 |  |
| c3.8 NSum machine | correct | 101 | correct | 3.91 | Mathilda 26x |
| c3.8 NSum wp40 | correct | 3497 | correct | 441 | Mathilda 7.9x |
| c3.9 NSum machine | correct | 6.6 | correct | 0.887 | Mathilda 7.4x |
| c3.9 NSum wp40 | correct | 24.2 | correct | 20.3 | Mathilda 1.2x |
| c3.10 NSum machine | correct | 14.7 | correct | 0.148 | Mathilda 1e+02x |
| c3.10 NSum wp40 | correct | 29.7 | correct | 12.5 | Mathilda 2.4x |
| c3.7 NSum machine | correct | 38.3 | correct | 2.34 | Mathilda 16x |
| c3.7 PolyLog wp40 | correct | 12.8 | correct | 0.502 | Mathilda 26x |
| c4.1 SeriesCoefficient exact | correct | 674 | correct | 102 | Mathilda 6.6x |
| c4.1 contour machine | correct | 1664 | correct | 44.6 | Mathilda 37x |
| c4.1 contour Trapezoidal wp40 | correct | 40410 | correct | 115377 | baseline 2.9x |
| c4.2 SeriesCoefficient exact | correct | 11.9 | correct | 4.25 | Mathilda 2.8x |
| c4.2 contour machine n=150 | correct | 305 | correct | 137 | Mathilda 2.2x |
| c4.2 contour Trapezoidal wp40 | correct | 1018 | none | - |  |
| c4.3 SeriesCoefficient exact | correct | 8.17 | correct | 4.73 | Mathilda 1.7x |
| c4.3 contour machine n=150 | correct | 155 | correct | 86.5 | Mathilda 1.8x |
| c4.3 contour Trapezoidal wp40 n=500 | correct | 811 | none | - |  |
| c4.4 Residue exact | correct | 111 | correct | 1.64 | Mathilda 68x |
| c4.4 contour machine | correct | 4.57 | wrong | 26 |  |
| c4.4 contour Trapezoidal wp40 | correct | 86.2 | correct | 394 | baseline 4.6x |
| c4.5 Residue exact | correct | 6.26 | correct | 3.92 | Mathilda 1.6x |
| c4.5 contour machine | correct | 6.96 | correct | 0.521 | Mathilda 13x |
| c4.5 contour Trapezoidal wp40 | correct | 44 | correct | 12.6 | Mathilda 3.5x |
| c4.6 SeriesCoefficient exact | correct | 5.12 | correct | 0.439 | Mathilda 12x |
| c4.6 contour machine | wrong | 2789 | wrong | 373 |  |
| c4.6 contour Trapezoidal wp40 | correct | 60270 | correct | 88733 | baseline 1.5x |
| c4.7 SeriesCoefficient exact | correct | 265 | correct | 24.8 | Mathilda 11x |
| c4.7 contour machine n=60 | correct | 19.6 | correct | 26.8 | baseline 1.4x |
| c4.7 contour Trapezoidal wp40 | correct | 184 | correct | 585 | baseline 3.2x |
| c4.8 Residue exact | correct | 7.66 | correct | 2.53 | Mathilda 3x |
| c4.8 contour machine | correct | 14.5 | correct | 10.5 | Mathilda 1.4x |
| c4.8 contour Trapezoidal wp40 | correct | 23.5 | correct | 54.9 | baseline 2.3x |
| c4.9 SeriesCoefficient exact | correct | 14.6 | correct | 141 | baseline 9.6x |
| c4.9 contour machine n=100 | correct | 10.2 | correct | 4.52 | Mathilda 2.3x |
| c4.9 contour Trapezoidal wp40 | correct | 24889 | correct | 127 | Mathilda 2e+02x |
| c4.10 Residue exact | correct | 4.47 | correct | 0.234 | Mathilda 19x |
| c4.10 contour machine | correct | 11.2 | correct | 5.64 | Mathilda 2x |
| c4.10 contour Trapezoidal wp40 | correct | 7.96 | correct | 18.3 | baseline 2.3x |
| c5.1 FindRoot machine | correct | 0.174 | correct | 0.032 | Mathilda 5.4x |
| c5.1 FindRoot wp40 | correct | 0.902 | correct | 0.306 | Mathilda 2.9x |
| c5.3 FindRoot machine | correct | 4.15 | correct | 0.067 | Mathilda 62x |
| c5.3 FindRoot wp40 | correct | 7.08 | correct | 0.538 | Mathilda 13x |
| c5.5 FindRoot machine | correct | 1.14 | none | - |  |
| c5.5 FindRoot wp40 | correct | 30 | none | - |  |
| c5.6 FindRoot machine | correct | 0.153 | correct | 0.021 | Mathilda 7.3x |
| c5.6 FindRoot wp40 | correct | 2.05 | correct | 3.7 | baseline 1.8x |
| c5.7 FindRoot machine | correct | 1.05 | correct | 0.028 | Mathilda 38x |
| c5.7 FindRoot wp40 | correct | 3.27 | correct | 0.786 | Mathilda 4.2x |
| c5.8 FindRoot machine | correct | 0.173 | correct | 0.033 | Mathilda 5.2x |
| c5.8 FindRoot wp40 | correct | 1.23 | correct | 0.434 | Mathilda 2.8x |
| c5.10 FindRoot machine | correct | 6.21 | none | - |  |
| c5.10 FindRoot wp40 | correct | 78.8 | none | - |  |
| c5.2 FindRoot machine | correct | 0.211 | correct | 0.084 | Mathilda 2.5x |
| c5.2 FindRoot wp40 | correct | 1.25 | correct | 0.543 | Mathilda 2.3x |
| c5.4 NSolve machine | wrong | 1.71 | correct | 0.414 |  |
| c5.4 NSolve wp60 | correct | 253 | correct | 0.599 | Mathilda 4.2e+02x |
| c5.9 NSolve machine | correct | 1.91 | correct | 1.93 | baseline 1x |
| c5.9 NSolve wp40 | correct | 15 | correct | 3.47 | Mathilda 4.3x |
| c6.1 GroebnerBasis exact | correct | 2.81 | correct | 0.811 | Mathilda 3.5x |
| c6.1 NSolve machine | correct | 21.8 | correct | 1.79 | Mathilda 12x |
| c6.1 NSolve wp40 | correct | 9.38 | correct | 1.26 | Mathilda 7.4x |
| c6.2 GroebnerBasis exact | correct | 0.947 | correct | 0.63 | Mathilda 1.5x |
| c6.2 NSolve machine | correct | 14.4 | correct | 0.452 | Mathilda 32x |
| c6.2 NSolve wp40 | correct | 12.6 | correct | 0.903 | Mathilda 14x |
| c6.3 GroebnerBasis exact | correct | 1.59 | correct | 0.762 | Mathilda 2.1x |
| c6.3 NSolve machine | correct | 16.6 | correct | 2.25 | Mathilda 7.4x |
| c6.3 NSolve wp40 | correct | 17.3 | correct | 3.19 | Mathilda 5.4x |
| c6.4 GroebnerBasis exact | correct | 1.52 | correct | 0.673 | Mathilda 2.3x |
| c6.4 NSolve machine | correct | 23.9 | correct | 2.49 | Mathilda 9.6x |
| c6.4 NSolve wp40 | correct | 31.1 | correct | 3.33 | Mathilda 9.4x |
| c6.5 GroebnerBasis exact | correct | 2.56 | correct | 0.846 | Mathilda 3x |
| c6.5 NSolve machine | correct | 26.7 | correct | 0.583 | Mathilda 46x |
| c6.5 NSolve wp40 | correct | 52.5 | correct | 2.06 | Mathilda 26x |
| c6.6 GroebnerBasis exact | correct | 15.6 | correct | 8.72 | Mathilda 1.8x |
| c6.6 NSolve machine | correct | 18.6 | correct | 4.38 | Mathilda 4.2x |
| c6.6 NSolve wp40 | correct | 447 | correct | 10.2 | Mathilda 44x |
| c6.7 GroebnerBasis exact | correct | 7.27 | correct | 0.922 | Mathilda 7.9x |
| c6.7 NSolve machine | correct | 73.6 | correct | 2.84 | Mathilda 26x |
| c6.7 NSolve wp40 | correct | 133 | correct | 5.41 | Mathilda 25x |
| c6.8 GroebnerBasis exact | correct | 0.808 | correct | 0.318 | Mathilda 2.5x |
| c6.8 NSolve machine | correct | 14.4 | correct | 0.63 | Mathilda 23x |
| c6.8 NSolve wp40 | correct | 18.1 | correct | 1.11 | Mathilda 16x |
| c6.9 GroebnerBasis exact | correct | 1.22 | correct | 0.616 | Mathilda 2x |
| c6.9 NSolve machine | correct | 16.5 | correct | 1.57 | Mathilda 10x |
| c6.9 NSolve wp40 | correct | 18.4 | correct | 2.44 | Mathilda 7.5x |
| c6.10 GroebnerBasis exact | correct | 9.48 | correct | 13.8 | baseline 1.5x |
| c6.10 NSolve machine | correct | 39.4 | correct | 9.88 | Mathilda 4x |
| c6.10 NSolve wp40 | correct | 46.4 | correct | 12.7 | Mathilda 3.7x |
| c7.1 Limit exact | correct | 272 | correct | 5071 | baseline 19x |
| c7.1 N machine | correct | 1.27 | correct | 0.117 | Mathilda 11x |
| c7.1 N wp30 | none | - | none | - |  |
| c7.2 Limit exact | correct | 29.6 | wrong | 5 |  |
| c7.2 N machine | correct | 0.094 | correct | 0.032 | Mathilda 2.9x |
| c7.2 N wp30 | wrong | 1.79 | correct | 0.921 |  |
| c7.3 Limit exact | correct | 22.6 | correct | 11.3 | Mathilda 2x |
| c7.3 N machine | correct | 0.079 | correct | 0.03 | Mathilda 2.6x |
| c7.3 N wp30 | correct | 0.621 | correct | 0.093 | Mathilda 6.7x |
| c7.4 Limit exact | correct | 14.5 | correct | 0.834 | Mathilda 17x |
| c7.4 N machine | correct | 0.054 | correct | 0.037 | Mathilda 1.5x |
| c7.4 N wp30 | correct | 0.281 | correct | 0.119 | Mathilda 2.4x |
| c7.5 Limit exact | correct | 27.9 | correct | 7.83 | Mathilda 3.6x |
| c7.5 N machine | correct | 0.07 | correct | 0.025 | Mathilda 2.8x |
| c7.5 N wp30 | correct | 0.624 | correct | 0.057 | Mathilda 11x |
| c7.6 Limit exact | correct | 58.9 | correct | 4.11 | Mathilda 14x |
| c7.6 N machine | correct | 0.074 | correct | 0.121 | baseline 1.6x |
| c7.6 N wp30 | correct | 0.328 | correct | 0.098 | Mathilda 3.3x |
| c7.7 Limit exact | correct | 10.8 | correct | 26.1 | baseline 2.4x |
| c7.7 N machine | correct | 0.056 | correct | 0.024 | Mathilda 2.3x |
| c7.7 N wp30 | correct | 0.163 | correct | 0.086 | Mathilda 1.9x |
| c7.8 Limit exact | correct | 277 | correct | 3.86 | Mathilda 72x |
| c7.8 N machine | correct | 0.027 | correct | 0.019 | Mathilda 1.4x |
| c7.8 N wp30 | correct | 1.43 | correct | 0.257 | Mathilda 5.6x |
| c7.9 Limit exact | correct | 34.9 | correct | 7.92 | Mathilda 4.4x |
| c7.9 N machine | correct | 0.084 | correct | 0.032 | Mathilda 2.6x |
| c7.9 N wp30 | correct | 0.198 | correct | 0.104 | Mathilda 1.9x |
| c7.10 Limit exact | correct | 318 | none | - |  |
| c7.10 N machine | correct | 0.109 | correct | 0.033 | Mathilda 3.3x |
| c7.10 N wp30 | correct | 1.03 | correct | 0.12 | Mathilda 8.6x |
| c8.1 Integrate exact | correct | 2267 | correct | 116 | Mathilda 20x |
| c8.1 NIntegrate machine | correct | 2.46 | correct | 0.054 | Mathilda 46x |
| c8.1 NIntegrate wp40 | correct | 24.6 | correct | 3.33 | Mathilda 7.4x |
| c8.2 Integrate exact | correct | 1100 | none | - |  |
| c8.2 NIntegrate machine | correct | 3.09 | correct | 0.046 | Mathilda 67x |
| c8.2 NIntegrate wp40 | correct | 58.3 | correct | 29.9 | Mathilda 1.9x |
| c8.3 Integrate exact | correct | 1137 | correct | 32.9 | Mathilda 35x |
| c8.3 NIntegrate machine | correct | 3.37 | correct | 0.03 | Mathilda 1.1e+02x |
| c8.3 NIntegrate wp40 | correct | 8.61 | correct | 3.76 | Mathilda 2.3x |
| c8.4 Integrate exact | none | - | none | - |  |
| c8.4 NIntegrate machine | correct | 3.26 | correct | 0.063 | Mathilda 52x |
| c8.4 NIntegrate wp40 | correct | 13.9 | correct | 8.62 | Mathilda 1.6x |
| c8.5 Integrate exact | correct | 2595 | correct | 2.71 | Mathilda 9.6e+02x |
| c8.5 NIntegrate machine | correct | 28 | correct | 24.3 | Mathilda 1.2x |
| c8.5 NIntegrate wp40 | correct | 1119 | correct | 849 | Mathilda 1.3x |
| c8.6 Integrate exact | none | - | none | - |  |
| c8.6 NIntegrate machine | correct | 2.94 | correct | 0.127 | Mathilda 23x |
| c8.6 NIntegrate wp40 | correct | 23.1 | correct | 4.35 | Mathilda 5.3x |
| c8.7 Integrate exact | correct | 442 | none | - |  |
| c8.7 NIntegrate machine | correct | 3.38 | correct | 0.028 | Mathilda 1.2e+02x |
| c8.7 NIntegrate wp40 | correct | 39.6 | correct | 62.3 | baseline 1.6x |
| c8.8 Integrate exact | correct | 11964 | correct | 37.7 | Mathilda 3.2e+02x |
| c8.8 NIntegrate machine | correct | 4.05 | correct | 0.055 | Mathilda 74x |
| c8.8 NIntegrate wp40 | correct | 15.7 | correct | 2.36 | Mathilda 6.7x |
| c8.9 Integrate exact | correct | 3143 | correct | 2.65 | Mathilda 1.2e+03x |
| c8.9 NIntegrate machine | correct | 4.15 | correct | 0.173 | Mathilda 24x |
| c8.9 NIntegrate wp40 | correct | 22.4 | correct | 2 | Mathilda 11x |
| c8.10 Integrate exact | correct | 1429 | correct | 115 | Mathilda 12x |
| c8.10 NIntegrate machine | correct | 2.93 | correct | 0.054 | Mathilda 54x |
| c8.10 NIntegrate wp40 | correct | 11.3 | correct | 4.23 | Mathilda 2.7x |
| c9.1 Eigenvalues machine n=6 | correct | 0.012 | correct | 0.016 | baseline 1.3x |
| c9.1 Eigenvalues wp300 n=60 | correct | 777 | correct | 524 | Mathilda 1.5x |
| c9.2 shooting machine | correct | 17.5 | correct | 2.36 | Mathilda 7.4x |
| c9.2 shooting wp40 | correct | 1442 | correct | 153 | Mathilda 9.4x |
| c9.3 Nystrom machine n=20 | correct | 2.79 | correct | 2.61 | Mathilda 1.1x |
| c9.3 Nystrom wp50 n=50 | correct | 1030 | correct | 112 | Mathilda 9.2x |
| c9.4 oscillator basis machine N=40 | correct | 9.81 | correct | 7.23 | Mathilda 1.4x |
| c9.4 oscillator basis wp60 N=100 | correct | 8877 | correct | 8314 | Mathilda 1.1x |
| c9.5 Eigenvalues machine | correct | 0.033 | correct | 0.013 | Mathilda 2.5x |
| c9.5 Eigenvalues wp40 | correct | 13.1 | correct | 0.71 | Mathilda 19x |
| c9.6 shooting machine | correct | 73 | correct | 2.93 | Mathilda 25x |
| c9.6 shooting wp40 | correct | 1395 | correct | 9.8 | Mathilda 1.4e+02x |
| c9.7 Eigenvalues machine N=20 | correct | 0.04 | correct | 0.054 | baseline 1.3x |
| c9.7 Eigenvalues wp60 N=40 | correct | 87.3 | correct | 9.71 | Mathilda 9x |
| c9.8 Eigenvalues machine n=12 | correct | 0.067 | correct | 0.051 | Mathilda 1.3x |
| c9.8 Eigenvalues wp200 n=40 | correct | 306 | correct | 119 | Mathilda 2.6x |
| c9.9 Eigenvalues machine n=8 | correct | 0.035 | correct | 0.038 | baseline 1.1x |
| c9.9 Eigenvalues wp200 n=30 | correct | 193 | correct | 61.8 | Mathilda 3.1x |
| c9.10 shooting machine | correct | 9.62 | correct | 0.814 | Mathilda 12x |
| c9.10 shooting wp40 | correct | 198 | correct | 8.32 | Mathilda 24x |
| c10.1 NDSolve pole machine | correct | 23.1 | correct | 1.46 | Mathilda 16x |
| c10.1 NDSolve pole wp40 | correct | 2406 | correct | 41.7 | Mathilda 58x |
| c10.2 NDSolve pole machine | correct | 3.95 | correct | 1.85 | Mathilda 2.1x |
| c10.2 NDSolve pole wp40 | correct | 2700 | correct | 55.5 | Mathilda 49x |
| c10.3 NDSolve pole machine | correct | 13.1 | correct | 1.77 | Mathilda 7.4x |
| c10.3 NDSolve pole wp40 | correct | 5189 | correct | 38.6 | Mathilda 1.3e+02x |
| c10.4 NDSolve stiff machine | correct | 5.47 | correct | 4.08 | Mathilda 1.3x |
| c10.4 NDSolve stiff wp32 | correct | 560 | correct | 26.8 | Mathilda 21x |
| c10.5 NDSolve BVP machine eps=1/100 | correct | 102 | none | - |  |
| c10.5 NDSolve BVP wp32 eps=1/1000 | correct | 2557 | none | - |  |
| c10.6 NDSolve chaotic machine T=50 | correct | 3.17 | correct | 5.04 | baseline 1.6x |
| c10.6 NDSolve chaotic wp40 T=100 | correct | 3018 | correct | 31 | Mathilda 97x |
| c10.7 NDSolve chaotic machine T=20 | correct | 3.22 | correct | 1.69 | Mathilda 1.9x |
| c10.7 NDSolve chaotic wp40 T=40 | correct | 1192 | correct | 30 | Mathilda 40x |
| c10.8 NDSolve singular start machine | correct | 1.82 | correct | 0.638 | Mathilda 2.9x |
| c10.8 NDSolve singular start wp40 | correct | 28.7 | correct | 1.7 | Mathilda 17x |
| c10.9 NDSolve machine | correct | 0.92 | correct | 0.172 | Mathilda 5.3x |
| c10.9 NDSolve wp40 | correct | 22 | correct | 2.52 | Mathilda 8.7x |
| c10.10 NDSolve backward machine L=8 | correct | 2.12 | none | - |  |
| c10.10 NDSolve backward wp40 L=20 | correct | 184 | none | - |  |
