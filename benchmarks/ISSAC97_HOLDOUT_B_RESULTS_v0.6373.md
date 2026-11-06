# ISSAC '97 holdout B (experiments 146-155): classification against the true answers

| family | cases | Mathilda correct | baseline correct | both correct | Mathilda faster | Mathilda slower | only Mathilda correct | only baseline correct | Mathilda wrong | Mathilda no answer |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 ill-conditioned linear algebra | 24 | 24 | 24 | 24 | 24 | 0 | 0 | 0 | 0 | 0 |
| 2 hard definite integrals | 21 | 21 | 16 | 16 | 16 | 0 | 5 | 0 | 0 | 0 |
| 3 slowly convergent series | 20 | 20 | 16 | 16 | 16 | 0 | 4 | 0 | 0 | 0 |
| 4 exact coefficients and residues | 30 | 30 | 27 | 27 | 27 | 0 | 3 | 0 | 0 | 0 |
| 5 zeros of special functions | 25 | 25 | 25 | 25 | 22 | 3 | 0 | 0 | 0 | 0 |
| 6 polynomial systems | 30 | 30 | 28 | 28 | 28 | 0 | 2 | 0 | 0 | 0 |
| 7 cancellation limits | 32 | 32 | 32 | 32 | 30 | 2 | 0 | 0 | 0 | 0 |
| 8 closed-form definite integrals | 30 | 30 | 29 | 29 | 29 | 0 | 1 | 0 | 0 | 0 |
| 9 operator eigenvalues | 22 | 22 | 22 | 22 | 18 | 4 | 0 | 0 | 0 | 0 |
| 10 ODE singularities | 24 | 24 | 24 | 24 | 24 | 0 | 0 | 0 | 0 | 0 |
| **all** | **258** | **258** | **243** | **243** | **234** | **9** | **15** | **0** | **0** | **0** |

## Cases

| case | Mathilda | baseline | Mathilda ms | baseline ms | ratio | speed |
|---|---|---|---|---|---|---|
| b1.1 Eigenvalues machine | CORRECT | CORRECT | 0.1 | 0.2 | 0.36 | faster |
| b1.1 Eigenvalues wp200 | CORRECT | CORRECT | 36.8 | 554.3 | 0.07 | faster |
| b1.2 SingularValueList machine | CORRECT | CORRECT | 0.1 | 0.1 | 0.48 | faster |
| b1.2 SingularValueList wp150 | CORRECT | CORRECT | 30.3 | 86.0 | 0.35 | faster |
| b1.3 LinearSolve exact | CORRECT | CORRECT | 2.5 | 31.2 | 0.08 | faster |
| b1.3 LinearSolve wp100 | CORRECT | CORRECT | 3.2 | 20.5 | 0.16 | faster |
| b1.3 LinearSolve machine | CORRECT | CORRECT | 0.2 | 0.5 | 0.33 | faster |
| b1.4 SingularValueList machine | CORRECT | CORRECT | 0.2 | 0.8 | 0.28 | faster |
| b1.4 SingularValueList wp100 | CORRECT | CORRECT | 288.1 | 1143.6 | 0.25 | faster |
| b1.5 Eigenvalues machine | CORRECT | CORRECT | 0.1 | 0.3 | 0.40 | faster |
| b1.5 Eigenvalues wp150 | CORRECT | CORRECT | 69.4 | 386.7 | 0.18 | faster |
| b1.6 LeastSquares exact | CORRECT | CORRECT | 16.4 | 48.6 | 0.34 | faster |
| b1.6 LeastSquares wp60 | CORRECT | CORRECT | 44.7 | 144.9 | 0.31 | faster |
| b1.6 LeastSquares machine | CORRECT | CORRECT | 0.4 | 1.2 | 0.36 | faster |
| b1.7 Det exact | CORRECT | CORRECT | 6.9 | 16.7 | 0.41 | faster |
| b1.7 Det wp120 | CORRECT | CORRECT | 1.8 | 6.5 | 0.27 | faster |
| b1.7 Det machine | CORRECT | CORRECT | 0.1 | 0.3 | 0.34 | faster |
| b1.8 Eigenvalues machine | CORRECT | CORRECT | 0.7 | 1.7 | 0.39 | faster |
| b1.8 Eigenvalues wp120 | CORRECT | CORRECT | 16.8 | 172.8 | 0.10 | faster |
| b1.9 Eigenvalues machine | CORRECT | CORRECT | 0.1 | 0.5 | 0.27 | faster |
| b1.9 Eigenvalues wp60 | CORRECT | CORRECT | 6.8 | 130.8 | 0.05 | faster |
| b1.10 Det exact | CORRECT | CORRECT | 4.0 | 9.7 | 0.41 | faster |
| b1.10 Det wp300 | CORRECT | CORRECT | 2.8 | 11.8 | 0.23 | faster |
| b1.10 Det machine | CORRECT | CORRECT | 0.0 | 0.1 | 0.46 | faster |
| b2.1 NIntegrate machine | CORRECT | CORRECT | 0.1 | 44.2 | 0.00 | faster |
| b2.1 NIntegrate wp30 | CORRECT | CORRECT | 4.4 | 715.8 | 0.01 | faster |
| b2.2 NIntegrate machine | CORRECT | WRONG | 8.3 |  |  |  |
| b2.2 NIntegrate wp30 | CORRECT | WRONG | 10.9 |  |  |  |
| b2.3 NIntegrate machine | CORRECT | CORRECT | 0.1 | 17.5 | 0.00 | faster |
| b2.3 NIntegrate wp30 | CORRECT | CORRECT | 2.1 | 153.9 | 0.01 | faster |
| b2.4 NIntegrate machine | CORRECT | CORRECT | 0.1 | 13.1 | 0.01 | faster |
| b2.4 NIntegrate wp30 | CORRECT | CORRECT | 4.6 | 48.7 | 0.09 | faster |
| b2.6 NIntegrate machine | CORRECT | CORRECT | 0.0 | 14.3 | 0.00 | faster |
| b2.6 NIntegrate wp30 | CORRECT | CORRECT | 1.8 | 22.5 | 0.08 | faster |
| b2.8 NIntegrate machine | CORRECT | CORRECT | 0.1 | 221.9 | 0.00 | faster |
| b2.8 NIntegrate wp30 | CORRECT | CORRECT | 5.4 | 588.1 | 0.01 | faster |
| b2.9 NIntegrate machine | CORRECT | WRONG | 6.3 |  |  |  |
| b2.9 NIntegrate wp30 | CORRECT | WRONG | 358.5 |  |  |  |
| b2.10 NIntegrate machine | CORRECT | CORRECT | 23.2 | 406.9 | 0.06 | faster |
| b2.10 NIntegrate wp30 | CORRECT | CORRECT | 579.5 | 10961.0 | 0.05 | faster |
| b2.5 NIntegrate machine | CORRECT | CORRECT | 0.9 | 39.0 | 0.02 | faster |
| b2.5 NIntegrate machine X=25 | CORRECT | CORRECT | 0.0 | 20.7 | 0.00 | faster |
| b2.5 NIntegrate wp30 | CORRECT | CORRECT | 3.1 | 69.1 | 0.05 | faster |
| b2.7 NIntegrate machine | CORRECT | CORRECT | 66.5 | 1219.0 | 0.05 | faster |
| b2.7 NIntegrate wp20 | CORRECT | NO ANSWER | 3050.0 |  |  |  |
| b3.1 NSum machine | CORRECT | CORRECT | 0.1 | 91.9 | 0.00 | faster |
| b3.1 NSum wp30 | CORRECT | CORRECT | 0.6 | 123.5 | 0.00 | faster |
| b3.2 NSum machine | CORRECT | CORRECT | 0.8 | 27.7 | 0.03 | faster |
| b3.2 NSum wp30 | CORRECT | CORRECT | 34.9 | 101.0 | 0.35 | faster |
| b3.3 NSum machine | CORRECT | CORRECT | 2.2 | 58.4 | 0.04 | faster |
| b3.3 NSum wp30 | CORRECT | CORRECT | 23.1 | 99.8 | 0.23 | faster |
| b3.4 NSum machine | CORRECT | CORRECT | 0.1 | 33.6 | 0.00 | faster |
| b3.4 NSum wp30 | CORRECT | CORRECT | 5.8 | 65.8 | 0.09 | faster |
| b3.5 NSum machine | CORRECT | CORRECT | 1.3 | 246.8 | 0.01 | faster |
| b3.5 NSum wp30 | CORRECT | CORRECT | 50.1 | 632.3 | 0.08 | faster |
| b3.6 NSum machine | CORRECT | CORRECT | 0.2 | 89.8 | 0.00 | faster |
| b3.6 NSum wp30 | CORRECT | CORRECT | 0.3 | 88.1 | 0.00 | faster |
| b3.7 NSum machine | CORRECT | CORRECT | 0.2 | 99.9 | 0.00 | faster |
| b3.7 NSum wp30 | CORRECT | CORRECT | 0.8 | 105.5 | 0.01 | faster |
| b3.8 NSum machine | CORRECT | WRONG | 0.0 |  |  |  |
| b3.8 NSum wp30 | CORRECT | WRONG | 0.0 |  |  |  |
| b3.9 NSum machine | CORRECT | WRONG | 0.0 |  |  |  |
| b3.9 NSum wp30 | CORRECT | WRONG | 0.0 |  |  |  |
| b3.10 NSum machine | CORRECT | CORRECT | 6.8 | 37904.7 | 0.00 | faster |
| b3.10 NSum wp30 | CORRECT | CORRECT | 1221.2 | 43368.9 | 0.03 | faster |
| b4.1 SeriesCoefficient exact | CORRECT | CORRECT | 899.4 | 34800.2 | 0.03 | faster |
| b4.1 contour wp60 | CORRECT | NO ANSWER | 10981.6 |  |  |  |
| b4.1 contour machine n=100 | CORRECT | CORRECT | 62.7 | 4654.1 | 0.01 | faster |
| b4.2 SeriesCoefficient exact | CORRECT | CORRECT | 61.8 | 3313.5 | 0.02 | faster |
| b4.2 contour wp160 | CORRECT | CORRECT | 54.8 | 649.3 | 0.08 | faster |
| b4.2 contour machine n=20 | CORRECT | CORRECT | 6.0 | 143.7 | 0.04 | faster |
| b4.3 Residue exact | CORRECT | CORRECT | 12.1 | 1168.9 | 0.01 | faster |
| b4.3 contour wp40 | CORRECT | CORRECT | 391.9 | 747.0 | 0.52 | faster |
| b4.3 contour machine | CORRECT | CORRECT | 2.1 | 344.4 | 0.01 | faster |
| b4.4 SeriesCoefficient exact | CORRECT | CORRECT | 39.9 | 3098.8 | 0.01 | faster |
| b4.4 contour wp80 | CORRECT | CORRECT | 10919.9 | 22417.0 | 0.49 | faster |
| b4.4 contour machine n=100 | CORRECT | CORRECT | 117.6 | 560.8 | 0.21 | faster |
| b4.5 SeriesCoefficient exact | CORRECT | CORRECT | 3.3 | 85.5 | 0.04 | faster |
| b4.5 contour wp100 | CORRECT | CORRECT | 490.1 | 767.9 | 0.64 | faster |
| b4.5 contour machine n=30 | CORRECT | CORRECT | 36.2 | 120.9 | 0.30 | faster |
| b4.6 SeriesCoefficient exact | CORRECT | CORRECT | 802.9 | 36328.8 | 0.02 | faster |
| b4.6 contour wp40 | CORRECT | NO ANSWER | 80949.9 |  |  |  |
| b4.6 contour machine | CORRECT | CORRECT | 6996.9 | 57545.8 | 0.12 | faster |
| b4.7 Residue exact | CORRECT | NO ANSWER | 0.1 |  |  |  |
| b4.7 contour wp50 | CORRECT | CORRECT | 12.0 | 62.6 | 0.19 | faster |
| b4.7 contour machine z^3 | CORRECT | CORRECT | 2.2 | 44.0 | 0.05 | faster |
| b4.8 Coefficient exact | CORRECT | CORRECT | 2.1 | 129.0 | 0.02 | faster |
| b4.8 contour wp30 | CORRECT | CORRECT | 45.7 | 318.4 | 0.14 | faster |
| b4.8 contour machine n=40 | CORRECT | CORRECT | 2.4 | 78.3 | 0.03 | faster |
| b4.9 SeriesCoefficient exact | CORRECT | CORRECT | 698.2 | 11677.9 | 0.06 | faster |
| b4.9 contour wp40 | CORRECT | CORRECT | 207.5 | 848.5 | 0.24 | faster |
| b4.9 contour machine n=30 | CORRECT | CORRECT | 6.4 | 295.4 | 0.02 | faster |
| b4.10 Residue exact | CORRECT | CORRECT | 127.1 | 160.1 | 0.79 | faster |
| b4.10 contour wp40 | CORRECT | CORRECT | 29.6 | 97.0 | 0.30 | faster |
| b4.10 contour machine z^10 | CORRECT | CORRECT | 5.9 | 121.9 | 0.05 | faster |
| b5.1 FindRoot wp60 | CORRECT | CORRECT | 111.2 | 122.5 | 0.91 | faster |
| b5.1 FindRoot machine n=100 | CORRECT | CORRECT | 0.1 | 1.3 | 0.04 | faster |
| b5.2 FindRoot wp40 | CORRECT | CORRECT | 2.3 | 46.3 | 0.05 | faster |
| b5.2 FindRoot machine | CORRECT | CORRECT | 0.1 | 0.4 | 0.15 | faster |
| b5.3 FindRoot wp40 | CORRECT | CORRECT | 0.2 | 0.8 | 0.21 | faster |
| b5.3 FindRoot machine | CORRECT | CORRECT | 0.0 | 0.3 | 0.07 | faster |
| b5.5 FindRoot wp40 | CORRECT | CORRECT | 0.2 | 0.8 | 0.22 | faster |
| b5.5 FindRoot machine | CORRECT | CORRECT | 0.0 | 0.2 | 0.06 | faster |
| b5.6 FindRoot wp40 | CORRECT | CORRECT | 4.4 | 1.9 | 2.34 | slower |
| b5.6 FindRoot machine | CORRECT | CORRECT | 0.1 | 0.5 | 0.30 | faster |
| b5.8 FindRoot wp40 | CORRECT | CORRECT | 8.0 | 47.5 | 0.17 | faster |
| b5.8 FindRoot machine | CORRECT | CORRECT | 0.1 | 8.0 | 0.01 | faster |
| b5.3 BesselJZero wp30 | CORRECT | CORRECT | 0.2 | 0.0 | 22.70 | slower |
| b5.4 FindRoot wp30 | CORRECT | CORRECT | 8.4 | 143.5 | 0.06 | faster |
| b5.4 FindRoot machine | CORRECT | CORRECT | 2.8 | 1.2 | 2.37 | slower |
| b5.4 ZetaZero wp30 | CORRECT | CORRECT | 6.6 | 148.4 | 0.04 | faster |
| b5.5 AiryAiZero wp30 | CORRECT | CORRECT | 0.1 | 2.3 | 0.03 | faster |
| b5.7 CountRoots exact | CORRECT | CORRECT | 0.3 | 3.0 | 0.10 | faster |
| b5.7 NSolve wp40 | CORRECT | CORRECT | 1.4 | 7.5 | 0.18 | faster |
| b5.7 FindRoot machine | CORRECT | CORRECT | 0.5 | 0.9 | 0.51 | faster |
| b5.9 CountRoots exact | CORRECT | CORRECT | 0.1 | 0.7 | 0.09 | faster |
| b5.9 NSolve wp80 | CORRECT | CORRECT | 4.6 | 6.7 | 0.68 | faster |
| b5.9 NSolve machine n=8 | CORRECT | CORRECT | 0.1 | 2.3 | 0.05 | faster |
| b5.10 FindRoot wp30 | CORRECT | CORRECT | 1.1 | 8.7 | 0.12 | faster |
| b5.10 FindRoot machine | CORRECT | CORRECT | 0.0 | 3.0 | 0.01 | faster |
| b6.1 GroebnerBasis exact | CORRECT | CORRECT | 14.4 | 1558.7 | 0.01 | faster |
| b6.1 NSolve wp40 | CORRECT | CORRECT | 25.5 | 3807.9 | 0.01 | faster |
| b6.1 NSolve machine | CORRECT | WRONG | 11.3 |  |  |  |
| b6.2 GroebnerBasis exact | CORRECT | CORRECT | 3.4 | 176.8 | 0.02 | faster |
| b6.2 NSolve wp40 | CORRECT | CORRECT | 12.0 | 1857.4 | 0.01 | faster |
| b6.2 NSolve machine | CORRECT | CORRECT | 5.8 | 995.8 | 0.01 | faster |
| b6.3 GroebnerBasis exact | CORRECT | CORRECT | 12.2 | 53.1 | 0.23 | faster |
| b6.3 NSolve wp40 | CORRECT | CORRECT | 10.2 | 425.6 | 0.02 | faster |
| b6.3 NSolve machine | CORRECT | CORRECT | 7.8 | 299.2 | 0.03 | faster |
| b6.4 GroebnerBasis exact | CORRECT | CORRECT | 11.2 | 73.6 | 0.15 | faster |
| b6.4 NSolve wp40 | CORRECT | CORRECT | 11.6 | 1144.3 | 0.01 | faster |
| b6.4 NSolve machine | CORRECT | CORRECT | 5.0 | 76.2 | 0.06 | faster |
| b6.7 GroebnerBasis exact | CORRECT | CORRECT | 75.0 | 719.3 | 0.10 | faster |
| b6.7 NSolve wp40 | CORRECT | CORRECT | 104.4 | 3752.4 | 0.03 | faster |
| b6.7 NSolve machine | CORRECT | CORRECT | 109.2 | 2064.9 | 0.05 | faster |
| b6.8 GroebnerBasis exact | CORRECT | CORRECT | 0.3 | 1.6 | 0.16 | faster |
| b6.8 NSolve wp40 | CORRECT | CORRECT | 5.4 | 360.8 | 0.02 | faster |
| b6.8 NSolve machine | CORRECT | CORRECT | 1.1 | 31.7 | 0.03 | faster |
| b6.9 GroebnerBasis exact | CORRECT | CORRECT | 0.2 | 1.9 | 0.10 | faster |
| b6.9 NSolve wp40 | CORRECT | CORRECT | 42.1 | 149.9 | 0.28 | faster |
| b6.9 NSolve machine | CORRECT | CORRECT | 15.3 | 83.3 | 0.18 | faster |
| b6.10 GroebnerBasis exact | CORRECT | CORRECT | 17.3 | 30.3 | 0.57 | faster |
| b6.10 NSolve wp40 | CORRECT | CORRECT | 14.5 | 162.0 | 0.09 | faster |
| b6.10 NSolve machine | CORRECT | WRONG | 10.8 |  |  |  |
| b6.5 GroebnerBasis exact | CORRECT | CORRECT | 0.3 | 1.9 | 0.14 | faster |
| b6.5 NSolve wp40 | CORRECT | CORRECT | 0.5 | 20.6 | 0.02 | faster |
| b6.5 NSolve machine | CORRECT | CORRECT | 0.4 | 10.4 | 0.04 | faster |
| b6.6 GroebnerBasis eliminate exact | CORRECT | CORRECT | 0.3 | 2.0 | 0.13 | faster |
| b6.6 NSolve wp40 | CORRECT | CORRECT | 2.6 | 125.6 | 0.02 | faster |
| b6.6 NSolve machine | CORRECT | CORRECT | 1.7 | 101.7 | 0.02 | faster |
| b7.1 Limit exact | CORRECT | CORRECT | 12.9 | 242.8 | 0.05 | faster |
| b7.1 approach wp300 | CORRECT | CORRECT | 3.4 | 3.8 | 0.91 | faster |
| b7.1 Neville machine | CORRECT | CORRECT | 0.2 | 0.6 | 0.27 | faster |
| b7.2 Limit exact | CORRECT | CORRECT | 74.2 | 191.7 | 0.39 | faster |
| b7.2 Neville wp600 | CORRECT | CORRECT | 1.9 | 1.7 | 1.08 | slower |
| b7.2 Neville machine | CORRECT | CORRECT | 0.0 | 0.1 | 0.52 | faster |
| b7.3 Limit exact | CORRECT | CORRECT | 1.3 | 208.4 | 0.01 | faster |
| b7.3 approach wp200 | CORRECT | CORRECT | 0.3 | 0.5 | 0.47 | faster |
| b7.3 Neville machine | CORRECT | CORRECT | 0.1 | 0.1 | 0.50 | faster |
| b7.4 Limit exact | CORRECT | CORRECT | 3.0 | 50.5 | 0.06 | faster |
| b7.4 approach wp100 | CORRECT | CORRECT | 1.2 | 5.3 | 0.23 | faster |
| b7.4 Neville machine | CORRECT | CORRECT | 0.1 | 0.3 | 0.26 | faster |
| b7.5 Limit exact | CORRECT | CORRECT | 6.3 | 91.2 | 0.07 | faster |
| b7.5 approach wp150 | CORRECT | CORRECT | 0.4 | 0.6 | 0.70 | faster |
| b7.5 Neville machine | CORRECT | CORRECT | 0.1 | 0.2 | 0.51 | faster |
| b7.6 Limit exact | CORRECT | CORRECT | 8.2 | 67.7 | 0.12 | faster |
| b7.6 approach wp150 | CORRECT | CORRECT | 0.6 | 0.9 | 0.61 | faster |
| b7.6 sample machine x=1e-8 | CORRECT | CORRECT | 0.0 | 0.0 | 0.24 | faster |
| b7.7 Limit exact | CORRECT | CORRECT | 4.2 | 102.7 | 0.04 | faster |
| b7.7 approach wp150 | CORRECT | CORRECT | 0.5 | 0.8 | 0.67 | faster |
| b7.7 Neville machine | CORRECT | CORRECT | 0.1 | 0.2 | 0.36 | faster |
| b7.8 Limit exact | CORRECT | CORRECT | 0.6 | 77.3 | 0.01 | faster |
| b7.8 approach wp250 | CORRECT | CORRECT | 0.3 | 1.2 | 0.23 | faster |
| b7.8 Neville machine | CORRECT | CORRECT | 0.1 | 0.2 | 0.22 | faster |
| b7.9 Limit exact | CORRECT | CORRECT | 1.2 | 93.8 | 0.01 | faster |
| b7.9 Neville wp100 | CORRECT | CORRECT | 2.8 | 1.7 | 1.63 | slower |
| b7.9 Neville machine | CORRECT | CORRECT | 0.1 | 0.2 | 0.50 | faster |
| b7.10 Limit exact | CORRECT | CORRECT | 97.9 | 164.0 | 0.60 | faster |
| b7.10 approach wp100 | CORRECT | CORRECT | 0.6 | 0.8 | 0.74 | faster |
| b7.10 Neville machine | CORRECT | CORRECT | 0.1 | 0.3 | 0.36 | faster |
| b7.3 Series exact | CORRECT | CORRECT | 0.6 | 17.7 | 0.04 | faster |
| b7.8 Series exact | CORRECT | CORRECT | 0.3 | 7.4 | 0.04 | faster |
| b8.1 Integrate exact | CORRECT | CORRECT | 66.5 | 2242.2 | 0.03 | faster |
| b8.1 NIntegrate wp40 | CORRECT | CORRECT | 2.7 | 13.7 | 0.20 | faster |
| b8.1 NIntegrate machine | CORRECT | CORRECT | 0.0 | 7.4 | 0.00 | faster |
| b8.2 Integrate exact | CORRECT | CORRECT | 82.8 | 4645.6 | 0.02 | faster |
| b8.2 NIntegrate wp40 | CORRECT | CORRECT | 3.6 | 54.1 | 0.07 | faster |
| b8.2 NIntegrate machine | CORRECT | CORRECT | 0.1 | 9.8 | 0.01 | faster |
| b8.3 Integrate exact | CORRECT | CORRECT | 26.0 | 1841.2 | 0.01 | faster |
| b8.3 NIntegrate wp40 | CORRECT | CORRECT | 2.4 | 11.2 | 0.21 | faster |
| b8.3 NIntegrate machine | CORRECT | CORRECT | 0.0 | 4.9 | 0.00 | faster |
| b8.4 Integrate exact | CORRECT | CORRECT | 5.6 | 3611.3 | 0.00 | faster |
| b8.4 NIntegrate wp40 | CORRECT | CORRECT | 14.2 | 60.6 | 0.23 | faster |
| b8.4 NIntegrate machine | CORRECT | CORRECT | 0.0 | 10.0 | 0.00 | faster |
| b8.5 Integrate exact | CORRECT | CORRECT | 6.8 | 14974.5 | 0.00 | faster |
| b8.5 NIntegrate wp40 | CORRECT | CORRECT | 10.7 | 1732.6 | 0.01 | faster |
| b8.5 NIntegrate machine | CORRECT | CORRECT | 0.1 | 30.7 | 0.00 | faster |
| b8.6 Integrate exact | CORRECT | CORRECT | 26.2 | 4078.2 | 0.01 | faster |
| b8.6 NIntegrate wp40 | CORRECT | CORRECT | 6.3 | 18.8 | 0.34 | faster |
| b8.6 NIntegrate machine | CORRECT | CORRECT | 0.0 | 7.1 | 0.01 | faster |
| b8.7 Integrate exact | CORRECT | CORRECT | 54.8 | 2409.3 | 0.02 | faster |
| b8.7 NIntegrate wp40 | CORRECT | CORRECT | 3.7 | 39.2 | 0.10 | faster |
| b8.7 NIntegrate machine | CORRECT | CORRECT | 0.1 | 8.7 | 0.01 | faster |
| b8.8 Integrate exact | CORRECT | CORRECT | 133.0 | 8299.6 | 0.02 | faster |
| b8.8 NIntegrate wp40 | CORRECT | CORRECT | 5.0 | 13.0 | 0.39 | faster |
| b8.8 NIntegrate machine | CORRECT | CORRECT | 0.0 | 6.5 | 0.00 | faster |
| b8.9 Integrate exact | CORRECT | CORRECT | 0.0 | 2643.4 | 0.00 | faster |
| b8.9 NIntegrate wp40 | CORRECT | CORRECT | 13.2 | 42.9 | 0.31 | faster |
| b8.9 NIntegrate machine | CORRECT | CORRECT | 0.0 | 8.2 | 0.00 | faster |
| b8.10 Integrate exact | CORRECT | NO ANSWER | 2.7 |  |  |  |
| b8.10 NIntegrate wp40 | CORRECT | CORRECT | 169.1 | 427.4 | 0.40 | faster |
| b8.10 NIntegrate machine | CORRECT | CORRECT | 0.1 | 18.9 | 0.00 | faster |
| b9.1 secular FindRoot wp40 | CORRECT | CORRECT | 0.1 | 0.6 | 0.14 | faster |
| b9.1 Nystrom Richardson machine | CORRECT | CORRECT | 578.2 | 988.8 | 0.58 | faster |
| b9.2 Gauss-Hermite Nystrom wp30 | CORRECT | CORRECT | 78.9 | 377.5 | 0.21 | faster |
| b9.2 Gauss-Hermite Nystrom machine | CORRECT | CORRECT | 14.3 | 46.8 | 0.31 | faster |
| b9.3 trapezoid Nystrom wp30 N=256 | CORRECT | CORRECT | 1907.8 | 17698.6 | 0.11 | faster |
| b9.3 trapezoid Nystrom machine N=256 | CORRECT | CORRECT | 343.0 | 585.8 | 0.59 | faster |
| b9.4 Eigenvalues wp50 N=150 | CORRECT | CORRECT | 182.4 | 4658.2 | 0.04 | faster |
| b9.4 Eigenvalues machine N=500 | CORRECT | CORRECT | 79.7 | 13.9 | 5.75 | slower |
| b9.5 Eigenvalues wp60 n=60 | CORRECT | CORRECT | 52.5 | 429.7 | 0.12 | faster |
| b9.5 Eigenvalues machine n=20 | CORRECT | CORRECT | 0.4 | 0.5 | 0.73 | faster |
| b9.6 Fourier matrix wp40 | CORRECT | CORRECT | 4.6 | 77.3 | 0.06 | faster |
| b9.6 Fourier matrix machine | CORRECT | CORRECT | 0.9 | 1.3 | 0.67 | faster |
| b9.6 MathieuCharacteristicA wp30 | CORRECT | CORRECT | 0.1 | 0.2 | 0.36 | faster |
| b9.6 MathieuCharacteristicA machine | CORRECT | CORRECT | 0.0 | 0.0 | 1.10 | slower |
| b9.7 harmonic basis wp50 N=100 | CORRECT | CORRECT | 82.2 | 1626.9 | 0.05 | faster |
| b9.7 harmonic basis machine N=60 | CORRECT | CORRECT | 9.7 | 9.1 | 1.07 | slower |
| b9.8 Gauss-Legendre Nystrom wp40 | CORRECT | CORRECT | 30.3 | 137.4 | 0.22 | faster |
| b9.8 Gauss-Legendre Nystrom machine | CORRECT | CORRECT | 5.7 | 7.4 | 0.77 | faster |
| b9.9 Eigenvalues wp40 N=10 | CORRECT | CORRECT | 47.0 | 591.0 | 0.08 | faster |
| b9.9 Eigenvalues machine N=40 | CORRECT | CORRECT | 217.0 | 174.8 | 1.24 | slower |
| b9.10 Eigenvalues pencil wp40 N=60 | CORRECT | CORRECT | 23.7 | 712.0 | 0.03 | faster |
| b9.10 Eigenvalues pencil machine N=400 | CORRECT | CORRECT | 403.0 | 457.6 | 0.88 | faster |
| b10.1 naive domain end machine | CORRECT | CORRECT | 0.4 | 2.3 | 0.18 | faster |
| b10.1 reciprocal switch wp40 | CORRECT | CORRECT | 2.3 | 75.0 | 0.03 | faster |
| b10.1 reciprocal switch machine | CORRECT | CORRECT | 0.2 | 1.9 | 0.11 | faster |
| b10.2 EllipticK exact | CORRECT | CORRECT | 0.0 | 0.0 | 0.64 | faster |
| b10.2 NDSolve quarter period wp40 | CORRECT | CORRECT | 2.5 | 65.9 | 0.04 | faster |
| b10.2 NDSolve quarter period machine | CORRECT | CORRECT | 0.2 | 1.1 | 0.18 | faster |
| b10.3 NDSolve x(20) wp40 | CORRECT | CORRECT | 54.4 | 1802.2 | 0.03 | faster |
| b10.3 NDSolve x(5) machine | CORRECT | CORRECT | 0.8 | 1.9 | 0.45 | faster |
| b10.4 series start + FindRoot wp30 | CORRECT | CORRECT | 2.1 | 23.3 | 0.09 | faster |
| b10.4 series start + FindRoot machine | CORRECT | CORRECT | 0.3 | 1.3 | 0.20 | faster |
| b10.5 Toepfer wp30 | CORRECT | CORRECT | 0.7 | 31.2 | 0.02 | faster |
| b10.5 Toepfer machine | CORRECT | CORRECT | 0.3 | 1.0 | 0.30 | faster |
| b10.6 NDSolve y(8) wp50 | CORRECT | CORRECT | 0.8 | 1173.2 | 0.00 | faster |
| b10.6 NDSolve y(3) machine | CORRECT | CORRECT | 0.1 | 0.8 | 0.14 | faster |
| b10.7 naive domain end machine | CORRECT | CORRECT | 0.7 | 2.9 | 0.23 | faster |
| b10.7 energy tail wp30 | CORRECT | CORRECT | 1.3 | 67.2 | 0.02 | faster |
| b10.7 energy tail machine | CORRECT | CORRECT | 0.2 | 1.2 | 0.20 | faster |
| b10.8 bisection shooting wp30 | CORRECT | CORRECT | 809.1 | 8711.7 | 0.09 | faster |
| b10.8 bisection shooting machine | CORRECT | CORRECT | 20.0 | 109.5 | 0.18 | faster |
| b10.9 naive domain end machine | CORRECT | CORRECT | 0.4 | 2.1 | 0.19 | faster |
| b10.9 inverse function wp30 | CORRECT | CORRECT | 3.5 | 88.5 | 0.04 | faster |
| b10.9 inverse function machine | CORRECT | CORRECT | 0.4 | 1.8 | 0.25 | faster |
| b10.10 NDSolve e=0.99 wp30 | CORRECT | CORRECT | 28.8 | 484.2 | 0.06 | faster |
| b10.10 NDSolve e=0.9 machine | CORRECT | CORRECT | 0.5 | 3.1 | 0.18 | faster |

## Baseline column: cases not answered correctly

| case | baseline | baseline check | true check |
|---|---|---|---|
| b2.2 NIntegrate machine | WRONG | `{539221, 2}` | `{999760, 2}` |
| b2.2 NIntegrate wp30 | WRONG | `{999755557848, 2}` | `{999760310630, 2}` |
| b2.9 NIntegrate machine | WRONG | `{48082, -3}` | `{62513, -3}` |
| b2.9 NIntegrate wp30 | WRONG | `{480819728211, -3}` | `{625129234764, -3}` |
| b2.7 NIntegrate wp20 | NO ANSWER | `"timeout"` | `{156725, 1}` |
| b3.8 NSum machine | WRONG | `{248120, 0}` | `{270581, 0}` |
| b3.8 NSum wp30 | WRONG | `{248120186455, 0}` | `{270580808428, 0}` |
| b3.9 NSum machine | WRONG | `{134158, 0}` | `{136843, 0}` |
| b3.9 NSum wp30 | WRONG | `{134157618212, 0}` | `{136843277762, 0}` |
| b4.1 contour wp60 | NO ANSWER | `"timeout"` | `{240614678640326, 31}` |
| b4.6 contour wp40 | NO ANSWER | `"timeout"` | `{-303284129702400, 16}` |
| b4.7 Residue exact | NO ANSWER | `"unevaluated"` | `{204822323834388, -20}` |
| b6.1 NSolve machine | WRONG | `{30, {33333333, -1}}` | `{32, {33333333, -1}}` |
| b6.10 NSolve machine | WRONG | `{9, {-67290476, -2}}` | `{10, {-67290476, -2}}` |
| b8.10 Integrate exact | NO ANSWER | `"unevaluated"` | `{384373946213433, -1}` |
